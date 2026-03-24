import json
import logging
import sys
from pathlib import Path

from rdkit import Chem

logger = logging.getLogger(__name__)


def _expand_path(path_value):
    if path_value is None:
        return None
    return str(Path(path_value).expanduser().resolve())


class _CarsiDockServer:
    def __init__(self):
        self.initialized = False
        self.model = None
        self.ligand_dict = None
        self.pocket_dict = None
        self.pocket = None
        self.docking = None
        self.read_ligands = None
        self.rtms_model = None
        self.rtms_pocket = None
        self.scoring = None
        self.device = None
        self.lbfgsbsrv = None
        self.args = None

    def init(self, payload):
        import torch
        carsidock_repo = _expand_path(payload["carsidock_repo"])
        if carsidock_repo not in sys.path:
            sys.path.insert(0, carsidock_repo)

        from src.utils.utils import get_carsidock_model, get_abs_path
        from src.utils.docking_inference_utils import read_ligands, docking
        from src.utils.docking_utils import extract_carsidock_pocket, extract_pocket
        from RTMScore.utils import scoring, get_rtmscore_model
        from run_screening import get_heavy_atom_positions

        self.read_ligands = read_ligands
        self.docking = docking
        self.scoring = scoring

        pdb_file = _expand_path(payload["pdb_file"])
        reflig = _expand_path(payload["reflig"])
        ckpt_path = _expand_path(payload["ckpt_path"])
        rtms_ckpt_path = _expand_path(payload["rtms_ckpt_path"])

        cuda_convert = bool(payload.get("cuda_convert", False))
        num_threads = int(payload.get("num_threads", 1))
        num_conformer = int(payload.get("num_conformer", 5))

        self.args = {
            "num_threads": num_threads,
            "num_conformer": num_conformer,
        }

        self.device = torch.device("cuda")

        if cuda_convert:
            import pydock
            self.lbfgsbsrv = pydock.LBFGSBServer(num_threads, 0)
        else:
            self.lbfgsbsrv = None

        self.model, self.ligand_dict, self.pocket_dict = get_carsidock_model(ckpt_path, self.device)
        self.rtms_model = get_rtmscore_model(rtms_ckpt_path)

        positions = get_heavy_atom_positions(reflig)
        self.pocket, _ = extract_carsidock_pocket(pdb_file, reflig)
        self.rtms_pocket = extract_pocket(pdb_file, positions, distance=10, del_water=True)

        self.initialized = True
        return {"status": "ok"}

    def dock_batch(self, payload):
        import torch
        if not self.initialized:
            return {"status": "error", "error": "Server not initialized"}

        import tempfile
        import os

        smiles_list = payload.get("ligands_smiles")
        mol_block_list = payload.get("ligands_mol_block")

        if not smiles_list and not mol_block_list:
            return {"status": "error", "error": "No ligands_smiles or ligands_mol_block provided"}

        ligands = smiles_list or mol_block_list
        use_smiles = smiles_list is not None

        molecule_id = payload.get("molecule_id")
        all_poses = []
        raw_blocks = []

        for ligand in ligands:
            try:
                if use_smiles:
                    init_mol_list = self.read_ligands(smiles=[ligand])[0]
                else:
                    mol = Chem.MolFromMolBlock(ligand, removeHs=False)
                    init_mol_list = self.read_ligands(mol_list=[mol])[0]
                torch.cuda.empty_cache()
                with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tmp:
                    tmp_path = tmp.name
                try:
                    outputs = self.docking(
                        self.model,
                        self.pocket,
                        init_mol_list,
                        self.ligand_dict,
                        self.pocket_dict,
                        device=self.device,
                        output_path=tmp_path,
                        num_threads=self.args["num_threads"],
                        lbfgsbsrv=self.lbfgsbsrv,
                    )
                    if os.path.exists(tmp_path):
                        with open(tmp_path) as f:
                            raw_blocks.append(f.read())
                finally:
                    if os.path.exists(tmp_path):
                        os.unlink(tmp_path)
                all_poses.extend(m for m in outputs["mol_list"] if m is not None)
            except Exception as e:
                logger.exception("CarsiDock docking failed for %s", molecule_id)
                return {"status": "error", "error": str(e)}

        all_poses = [Chem.RemoveHs(m) for m in all_poses]
        if not all_poses:
            return {"status": "ok", "results": []}

        best_pose = all_poses[0]
        best_pose.SetProp('_Name', molecule_id)
        _, rtms_scores = self.scoring(self.rtms_pocket, [best_pose], self.rtms_model)
        raw_block = "".join(raw_blocks)

        return {"status": "ok", "results": [{
            "name": molecule_id,
            "docking_score": float(rtms_scores[0]),
            "mol_block": Chem.MolToMolBlock(best_pose),
            "raw_block": raw_block,
        }]}


def _emit_server_response(req_id, payload):
    response = {"id": req_id, "payload": payload}
    sys.stdout.write(json.dumps(response) + "\n")
    sys.stdout.flush()


def run_carsidock_server():
    logging.basicConfig(level=logging.INFO)
    server = _CarsiDockServer()

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue

        try:
            request = json.loads(line)
            req_id = request.get("id")
            command = request.get("command")
            payload = request.get("payload") or {}

            if command == "init":
                result = server.init(payload)
            elif command == "dock_batch":
                result = server.dock_batch(payload)
            else:
                result = {"status": "error", "error": f"Unknown command: {command}"}

        except Exception as e:
            result = {"status": "error", "error": str(e)}
            req_id = request.get("id") if "request" in locals() else None

        _emit_server_response(req_id, result)

    return 0


if __name__ == "__main__":
    raise SystemExit(run_carsidock_server())
