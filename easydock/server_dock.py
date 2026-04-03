import atexit
import logging
import os
import threading
import timeit
from pathlib import Path
from typing import List

import yaml
from rdkit import Chem

from easydock.containers import build_server_container_cmd
from easydock.persistent_client import JsonLineProcessClient

logger = logging.getLogger(__name__)

ligand_preparation = None
pdbqt2molblock = None

_worker_client = None
_worker_key = None
_worker_lock = threading.Lock()


def _safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _normalize_path_values(obj):
    if isinstance(obj, dict):
        return {k: _normalize_path_values(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_normalize_path_values(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    return obj


def _ensure_preparation_functions():
    global ligand_preparation, pdbqt2molblock
    if ligand_preparation is None or pdbqt2molblock is None:
        from easydock.preparation_for_docking import (
            ligand_preparation as _ligand_preparation,
            pdbqt2molblock as _pdbqt2molblock,
        )
        if ligand_preparation is None:
            ligand_preparation = _ligand_preparation
        if pdbqt2molblock is None:
            pdbqt2molblock = _pdbqt2molblock

    return ligand_preparation, pdbqt2molblock


def _collect_path_dirs(payload):
    """Recursively collect unique parent dirs of all existing file/dir paths in payload."""
    dirs = set()
    items = (payload.values() if isinstance(payload, dict) else
             payload if isinstance(payload, (list, tuple)) else [payload])
    for v in items:
        if isinstance(v, str):
            p = os.path.abspath(os.path.expanduser(v))
            if os.path.isfile(p):
                dirs.add(os.path.dirname(p))
            elif os.path.isdir(p):
                dirs.add(p)
        elif isinstance(v, (dict, list)):
            dirs.update(_collect_path_dirs(v))
    return sorted(dirs)


def _parse_config(config_fname):
    with open(config_fname) as f:
        config = yaml.safe_load(f) or {}

    if not isinstance(config, dict):
        raise ValueError("Docking server config must be a mapping")

    if "script_file" not in config:
        raise ValueError("'script_file' is required in server docking config")

    config["script_file"] = str(config["script_file"])
    if "server_cwd" in config and config["server_cwd"] is not None:
        config["server_cwd"] = str(config["server_cwd"])

    config.setdefault("score_mode", "min")
    config.setdefault("startup_timeout", 120)
    config.setdefault("request_timeout", None)
    config.setdefault("boron_replacement", False)
    config.setdefault("ligand_in_format", "pdbqt")
    config.setdefault("ligand_out_format", "pdbqt")

    control_keys = {
        "script_file",
        "server_cwd",
        "score_mode",
        "startup_timeout",
        "request_timeout",
        "boron_replacement",
        "ligand_in_format",
        "ligand_out_format",
        "init_server",
    }

    init_server = config.get("init_server")
    if not isinstance(init_server, dict):
        init_server = {
            k: _normalize_path_values(v)
            for k, v in config.items()
            if k not in control_keys
        }

    config["_resolved_init_server"] = init_server
    config["_worker_key"] = yaml.safe_dump(
        {
            "script_file": config["script_file"],
            "server_cwd": config.get("server_cwd"),
            "request_timeout": config.get("request_timeout"),
            "init_server": init_server,
        },
        sort_keys=True,
    )

    bind_dirs = _collect_path_dirs(init_server)
    config["_launch_command"] = build_server_container_cmd(config["script_file"], bind_dirs)

    return config


def _prepare_ligand_payload(mol, config, ring_sample=False):
    """
    Return a list of instances in a specific format. Commonly there will be a single instance per molecule.
    Instances can be conformers if ring_sample is True
    :param mol:
    :param config:
    :param ring_sample:
    :return:
    """


    payload_type = config["ligand_in_format"]

    if payload_type == "pdbqt":
        ligand_preparation, _ = _ensure_preparation_functions()
        ligand_payload = ligand_preparation(
            mol,
            boron_replacement=False,
            ring_sample=ring_sample,
        )
        return ligand_payload

    if payload_type == "smiles":
        return [Chem.MolToSmiles(mol, isomericSmiles=True)]

    if payload_type == "mol":
        return [Chem.MolToMolBlock(mol)]

    raise ValueError(f"Unsupported ligand_in_format: {payload_type}")


def _response_payload(response):
    if isinstance(response, dict):
        payload = response.get("payload")
        if isinstance(payload, dict):
            return payload
    return {}


def _check_response_ok(response, context="request"):
    payload = _response_payload(response)
    status = payload.get("status", response.get("status") if isinstance(response, dict) else None)
    if status in (None, True, "ok", "OK"):
        return

    error = payload.get("error") or (
        response.get("error") if isinstance(response, dict) else None
    )
    if error:
        raise RuntimeError(f"{context} failed: {error}")
    raise RuntimeError(f"{context} failed with status={status!r}")


def _create_client(config):
    client = JsonLineProcessClient(
        command=config["_launch_command"],
        startup_timeout=config.get("startup_timeout", 120),
        request_timeout=config.get("request_timeout"),
        cwd=config.get("server_cwd"),
    )

    response = client.request(
        command="init",
        payload=config["_resolved_init_server"],
        timeout=config.get("request_timeout"),
    )
    _check_response_ok(response, context="init")
    return client


def _get_worker_client(config):
    global _worker_client, _worker_key

    with _worker_lock:
        need_new = (
            _worker_client is None
            or not _worker_client.is_alive()
            or _worker_key != config["_worker_key"]
        )

        if need_new:
            if _worker_client is not None:
                try:
                    _worker_client.close()
                except Exception:
                    logger.exception("Failed to close previous worker client")
            _worker_client = _create_client(config)
            _worker_key = config["_worker_key"]

        return _worker_client


def _extract_batch_results(response):
    """Parse response into a list of (molecule_id, result_items) tuples.

    The outer response payload contains a list under "results". Each element
    of that list is a per-molecule dict with molecule_id, status, error, and its own
    "results" list of "conformer" dicts.
    Returns a list of (molecule_id, conformer_list) where conformer_list is None on error.
    """
    if not isinstance(response, dict):
        return []

    payload = _response_payload(response)
    per_mol_list = payload.get("results")

    if not isinstance(per_mol_list, list):
        logger.warning("Received response for 'payload' -> 'results' is not a list. This is a break of API.")
        return []

    outputs = []
    for item in per_mol_list:
        if not isinstance(item, dict):
            logger.warning("Per-molecule result item is not a dict: %r", item)
            continue

        mol_id = item.get("molecule_id")
        status = item.get("status")
        error = item.get("error")

        if error or status not in (None, True, "ok", "OK"):
            logger.warning("Docking failed for %s: %s", mol_id, error or f"status={status!r}")
            outputs.append((mol_id, None))
            continue

        conformers = item.get("results")
        if not isinstance(conformers, list):
            logger.warning("Missing or invalid 'results' in per-molecule result for %s", mol_id)
            outputs.append((mol_id, None))
            continue

        outputs.append((mol_id, conformers))

    return outputs


def _extract_molecule_name(response):
    payload = _response_payload(response)
    name = payload.get("molecule_id")
    if name is None and isinstance(response, dict):
        name = response.get("molecule_id")
    return name if isinstance(name, str) and name.strip() else None


def _set_mol_block_name(mol_block, name):
    first_newline = mol_block.find('\n')
    if first_newline == -1:
        return mol_block
    return name + mol_block[first_newline:]


def _choose_best(items, score_mode="min"):
    if not items:
        return None
    if score_mode == "max":
        return max(items, key=lambda x: x["docking_score"])
    return min(items, key=lambda x: x["docking_score"])


def mol_dock(mols: Chem.Mol | List[Chem.Mol], config, ring_sample=False):

    if isinstance(mols, Chem.Mol):
        mols = [mols]

    config = _parse_config(config)

    data = []
    for mol in mols:
        mol_id = mol.GetProp("_Name")
        try:
            ligand_payload = _prepare_ligand_payload(mol, config, ring_sample=ring_sample)
            data.append((mol_id, ligand_payload))
        except Exception as e:
            logger.warning("Ligand preparation failed for %s: %s", mol_id, e)
            data.append((mol_id, None))
    if all(payload is None for mol_id, payload in data):
        return data  # list of (mol_id, None)

    failed_mol_ids = [mol_id for mol_id, ligand_payload in data if ligand_payload is None]
    data = [item for item in data if item[1] is not None]  # keep only prepared mols

    start_time = timeit.default_timer()

    try:
        client = _get_worker_client(config)

        payload = []
        for mol_id, ligand_payload in data:
            payload.append({
                "molecule_id": mol_id,
                # "ring_sample": bool(ring_sample),
                config["ligand_in_format"]: ligand_payload,
            })

        response = client.request(
            command="dock",
            payload=payload,
            timeout=config.get("request_timeout"),
        )
        _check_response_ok(response, context="dock")
    except Exception as e:
        mol_ids = [mol_id for mol_id, _ in data]
        logger.warning("Docking request failed for the whole batch of molecules %s: %s",
                       ', '.join(mol_ids), e)
        return [(mol_id, None) for mol_id in mol_ids + failed_mol_ids]

    dock_time = round((timeit.default_timer() - start_time) / len(data), 1)  # average time per mol
    results = _extract_batch_results(response)

    score_key = "docking_score"
    pose_key = "raw_block"
    mol_block_key = "mol_block"
    out_format = config.get("ligand_out_format", "pdbqt")
    
    score_mode = config.get("score_mode", "min")
    mol_by_id = {mol.GetProp("_Name"): mol for mol in mols}

    outputs = [(mol_id, None) for mol_id in failed_mol_ids]

    for mol_id, result_list in results:
        
        if result_list is None:
            outputs.append((mol_id, None))
            continue

        mol = mol_by_id.get(mol_id)
        dock_output_conformer_list = []

        for item in result_list:
            docking_score = _safe_float(item.get(score_key))
            if docking_score is None:
                continue

            mol_block = item.get(mol_block_key)
            raw_block = item.get(pose_key)

            if isinstance(mol_block, str) and mol_block.strip():
                entry = {"docking_score": docking_score, "mol_block": mol_block}
                if isinstance(raw_block, str):
                    entry["raw_block"] = raw_block
                dock_output_conformer_list.append(entry)
                continue

            if not isinstance(raw_block, str) or not raw_block.strip():
                continue

            # parse raw_block if mol_block is absent
            if out_format == "pdbqt":
                if "MODEL" not in raw_block:
                    continue
                try:
                    _, _pdbqt2molblock = _ensure_preparation_functions()
                    mol_block = _pdbqt2molblock(raw_block.split("MODEL", 1)[1], mol, mol_id)
                except Exception:
                    logger.exception("Failed to convert pdbqt pose for %s", mol_id)
                    continue
                dock_output_conformer_list.append({"docking_score": docking_score, "mol_block": mol_block, "raw_block": raw_block})

            elif out_format == "sdf":
                dock_output_conformer_list.append({"docking_score": docking_score, "mol_block": raw_block.split('$$$$\n')[0], "raw_block": raw_block})

            else:
                raise NotImplementedError(f'Output raw docking format {out_format!r} is not implemented.')

        if dock_output_conformer_list:
            output = _choose_best(dock_output_conformer_list, score_mode)
            output["mol_block"] = _set_mol_block_name(output["mol_block"], mol_id)
            output["dock_time"] = dock_time
        else:
            output = None

        outputs.append((mol_id, output))

    return outputs


def _cleanup_worker_client():
    global _worker_client
    if _worker_client is not None:
        try:
            _worker_client.close()
        except Exception:
            pass


atexit.register(_cleanup_worker_client)


def _reset_worker_state_for_tests():
    global _worker_client, _worker_key
    with _worker_lock:
        if _worker_client is not None:
            try:
                _worker_client.close()
            except Exception:
                pass
        _worker_client = None
        _worker_key = None
