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
        self.num_threads = 1
        self.num_conformer = 5

    def init(self, payload):
        import torch

        sys.path.insert(0, '/carsidock')

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

        ckpt_path = '/carsidock/checkpoints/carsidock_230731.ckpt'
        rtms_ckpt_path = '/carsidock/checkpoints/rtmscore_model1.pth'

        cuda_convert = bool(payload.get("cuda_convert", True))
        self.num_threads = int(payload.get("num_threads", self.num_threads))
        self.num_conformer = int(payload.get("num_conformer", self.num_conformer))

        self.device =  torch.device(payload.get("device", "cuda"))

        if cuda_convert:
            import pydock
            self.lbfgsbsrv = pydock.LBFGSBServer(self.num_threads, 0)
        else:
            self.lbfgsbsrv = None

        self.model, self.ligand_dict, self.pocket_dict = get_carsidock_model(ckpt_path, self.device)
        self.rtms_model = get_rtmscore_model(rtms_ckpt_path)

        positions = get_heavy_atom_positions(reflig)
        self.pocket, _ = extract_carsidock_pocket(pdb_file, reflig)
        self.rtms_pocket = extract_pocket(pdb_file, positions, distance=10, del_water=True)

        self.initialized = True
        return {"status": "ok"}

    def dock(self, payload):
        import torch
        if not self.initialized:
            return {"status": "error", "error": "Server not initialized"}

        import tempfile
        import os

        if not isinstance(payload, list):
            return {"status": "error", "error": f"payload must be a list of molecule dicts. Payload: {payload}"}

        nconf = max(self.num_conformer, 5)
        per_mol_results = []

        for mol_entry in payload:
            if not isinstance(mol_entry, dict):
                logger.warning("Skipping non-dict entry in payload: %r", mol_entry)
                continue

            molecule_id = mol_entry.get("molecule_id")
            smiles_list = mol_entry.get("smiles")
            mol_block_list = mol_entry.get("mol_block")

            if not smiles_list and not mol_block_list:
                per_mol_results.append({
                    "molecule_id": molecule_id,
                    "status": "error",
                    "error": "No smiles or mol_block provided",
                    "results": [],
                })
                continue

            ligands = smiles_list if smiles_list is not None else mol_block_list
            use_smiles = smiles_list is not None
            mol_output = []

            try:
                for ligand in ligands:
                    raw_block = None
                    if use_smiles:
                        init_mol_list = self.read_ligands(smiles=[ligand], num_use_conf=nconf)[0]
                    else:
                        mol = Chem.MolFromMolBlock(ligand, removeHs=False)
                        init_mol_list = self.read_ligands(mol_list=[mol], num_use_conf=nconf)[0]
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
                            num_threads=self.num_threads,
                            lbfgsbsrv=self.lbfgsbsrv,
                        )
                        if os.path.exists(tmp_path):
                            with open(tmp_path) as f:
                                raw_block = f.read()
                                raw_block = '$$$$\n'.join(raw_block.split('$$$$\n', self.num_conformer)[:self.num_conformer])
                    finally:
                        if os.path.exists(tmp_path):
                            os.unlink(tmp_path)

                    best_pose = outputs["mol_list"][0]
                    _, rtms_scores = self.scoring(self.rtms_pocket, [best_pose], self.rtms_model)
                    mol_output.append({
                        "docking_score": float(rtms_scores[0]),
                        "mol_block": Chem.MolToMolBlock(best_pose),
                        "raw_block": raw_block,
                    })

                per_mol_results.append({
                    "molecule_id": molecule_id,
                    "status": "ok",
                    "results": mol_output,
                })

            except Exception as e:
                logger.exception("CarsiDock docking failed for %s", molecule_id)
                per_mol_results.append({
                    "molecule_id": molecule_id,
                    "status": "error",
                    "error": str(e),
                    "results": [],
                })

        return {"status": "ok", "results": per_mol_results}


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
            elif command == "dock":
                result = server.dock(payload)
            else:
                result = {"status": "error", "error": f"Unknown command: {command}"}

        except Exception as e:
            result = {"status": "error", "error": str(e)}
            req_id = request.get("id") if "request" in locals() else None

        _emit_server_response(req_id, result)

    return 0


if __name__ == "__main__":
    raise SystemExit(run_carsidock_server())
