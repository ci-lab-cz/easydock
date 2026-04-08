import json
import logging
import sys
from pathlib import Path

from rdkit import Chem

logger = logging.getLogger(__name__)

INFO = {
    'batch_size': 10,
    'ligand_in_format': 'smiles',
    'ligand_out_format': 'sdf',
    'score_mode': 'max',
}


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
        pdb_file = _expand_path(payload.get("protein"))
        reflig = _expand_path(payload.get("reflig"))

        if not pdb_file:
            return {"status": "error", "error": "'protein' is required in init_server"}
        if not reflig:
            return {"status": "error", "error": "'reflig' is required in init_server"}

        cuda_convert = bool(payload.get("cuda_convert", True))
        self.num_threads = int(payload.get("num_threads", self.num_threads))
        self.num_conformer = int(payload.get("num_conformer", self.num_conformer))

        try:
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

            self.device = torch.device(payload.get("device", "cuda"))

            if cuda_convert:
                import pydock
                self.lbfgsbsrv = pydock.LBFGSBServer(self.num_threads, 0)
            else:
                self.lbfgsbsrv = None

            ckpt_path = '/carsidock/checkpoints/carsidock_230731.ckpt'
            rtms_ckpt_path = '/carsidock/checkpoints/rtmscore_model1.pth'

            self.model, self.ligand_dict, self.pocket_dict = get_carsidock_model(ckpt_path, self.device)
            self.rtms_model = get_rtmscore_model(rtms_ckpt_path)

            positions = get_heavy_atom_positions(reflig)
            self.pocket, _ = extract_carsidock_pocket(pdb_file, reflig)
            self.rtms_pocket = extract_pocket(pdb_file, positions, distance=10, del_water=True)

        except Exception as e:
            logger.exception("CarsiDock initialization failed")
            return {"status": "error", "error": str(e)}

        self.initialized = True
        return {"status": "ok"}

    def dock(self, payload):
        import torch
        if not self.initialized:
            return {"status": "error", "error": "Server not initialized"}

        import tempfile
        import os

        if not isinstance(payload, dict):
            return {"status": "error", "error": f"payload must be a dict of {{molecule_name: representation}}. Got: {type(payload)}"}

        nconf = max(self.num_conformer, 5)
        results = {}

        for molecule_id, instance in payload.items():
            if not instance:
                logger.warning("Empty instance for %s", molecule_id)
                results[molecule_id] = None
                continue

            try:
                raw_block = None
                if INFO['ligand_in_format'] == 'smiles':
                    init_mol_list = self.read_ligands(smiles=[instance], num_use_conf=nconf)[0]
                elif INFO['ligand_in_format'] == 'mol':
                    mol = Chem.MolFromMolBlock(instance, removeHs=False)
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
                results[molecule_id] = {
                    "docking_score": float(rtms_scores[0]),
                    "mol_block": Chem.MolToMolBlock(best_pose),
                    "raw_block": raw_block,
                }

            except Exception as e:
                logger.exception("CarsiDock docking failed for %s", molecule_id)
                results[molecule_id] = None

        return {"status": "ok", "results": results}


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

            if command == "info":
                result = INFO.copy()
            elif command == "init":
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
