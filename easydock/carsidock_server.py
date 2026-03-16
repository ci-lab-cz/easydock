import json
import logging
import sys
from pathlib import Path

from rdkit import Chem

logger = logging.getLogger(__name__)


def _safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _expand_path(path_value):
    if path_value is None:
        return None
    return str(Path(path_value).expanduser().resolve())


def _score_from_mol(mol):
    if mol is None:
        return None

    for key in ("docking_score", "score", "affinity", "loss"):
        if mol.HasProp(key):
            score = _safe_float(mol.GetProp(key))
            if score is not None:
                return score
    return None


def _coerce_score_from_outputs(output_item, mol):
    if isinstance(output_item, dict):
        for key in ("docking_score", "score", "affinity", "loss"):
            score = _safe_float(output_item.get(key))
            if score is not None:
                return score
    return _score_from_mol(mol)


class _CarsiDockServer:
    def __init__(self):
        self.initialized = False
        self.model = None
        self.ligand_dict = None
        self.pocket_dict = None
        self.pocket = None
        self.docking = None
        self.read_ligands = None
        self.args = None

    def init(self, payload):
        carsidock_repo = _expand_path(payload["carsidock_repo"])
        if carsidock_repo not in sys.path:
            sys.path.insert(0, carsidock_repo)

        from run_screening import get_carsidock_model
        from src.utils.docking_inference_utils import (
            extract_carsidock_pocket,
            read_ligands,
            prepare_data_from_mol,
        )

        self.read_ligands = read_ligands
        self.docking = prepare_data_from_mol

        pdb_file = _expand_path(payload["pdb_file"])
        sdf_file = _expand_path(payload["sdf_file"])
        ckpt_path = _expand_path(payload["ckpt_path"])

        cuda_convert = bool(payload.get("cuda_convert", False))
        num_threads = int(payload.get("num_threads", 1))
        num_conformer = int(payload.get("num_conformer", 5))

        self.args = {
            "pdb_file": pdb_file,
            "sdf_file": sdf_file,
            "ckpt_path": ckpt_path,
            "cuda_convert": cuda_convert,
            "num_threads": num_threads,
            "num_conformer": num_conformer,
        }

        (
            self.model,
            self.ligand_dict,
            self.pocket_dict,
        ) = get_carsidock_model(ckpt_path, device="cuda")

        self.pocket = extract_carsidock_pocket(
            pdb_file,
            sdf_file,
            threading=num_threads,
            cuda_convert=cuda_convert,
        )

        self.initialized = True
        return {"status": "ok"}

    def dock_batch(self, payload):
        if not self.initialized:
            return {"status": "error", "error": "Server not initialized"}

        smiles_list = payload.get("ligands_smiles")
        if not smiles_list:
            return {"status": "error", "error": "No ligands_smiles provided"}

        molecule_id = payload.get("molecule_id")
        results = []

        for smiles in smiles_list:
            try:
                init_mol_list, _ = self.read_ligands(
                    smiles=[smiles],
                    num_use_conf=self.args["num_conformer"],
                )

                output_list, _ = self.docking(
                    self.model,
                    self.pocket,
                    init_mol_list,
                    self.ligand_dict,
                    self.pocket_dict,
                    device="cuda",
                    bsz=1,
                )

                best_item = None
                for output_item in output_list:
                    mol = None
                    if isinstance(output_item, dict):
                        mol = output_item.get("mol") or output_item.get("rdmol")
                    score = _coerce_score_from_outputs(output_item, mol)
                    if mol is None or score is None:
                        continue

                    mol_block = Chem.MolToMolBlock(mol)
                    candidate = {
                        "name": molecule_id,
                        "docking_score": score,
                        "mol_block": mol_block,
                    }
                    if best_item is None or score < best_item["docking_score"]:
                        best_item = candidate

                if best_item is not None:
                    results.append(best_item)

            except Exception as e:
                logger.exception("CarsiDock docking failed for %s", molecule_id)
                return {"status": "error", "error": str(e)}

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
