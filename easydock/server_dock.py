import atexit
import logging
import threading
import timeit
from pathlib import Path

import yaml
from rdkit import Chem

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

    config.setdefault("init_command", "init")
    config.setdefault("dock_command", "dock")
    config.setdefault("result_items_key", "results")
    config.setdefault("score_key", "docking_score")
    config.setdefault("pose_key", "raw_block")
    config.setdefault("mol_block_key", "mol_block")
    config.setdefault("score_mode", "min")
    config.setdefault("startup_timeout", 120)
    config.setdefault("request_timeout", None)
    config.setdefault("boron_replacement", False)
    config.setdefault("ligand_payload_type", "pdbqt")
    config.setdefault("raw_format", "pdbqt")

    control_keys = {
        "script_file",
        "server_cwd",
        "init_command",
        "dock_command",
        "result_items_key",
        "score_key",
        "pose_key",
        "mol_block_key",
        "score_mode",
        "startup_timeout",
        "request_timeout",
        "boron_replacement",
        "ligand_payload_type",
        "raw_format",
        "init_payload",
    }

    init_payload = config.get("init_payload")
    if not isinstance(init_payload, dict):
        init_payload = {
            k: _normalize_path_values(v)
            for k, v in config.items()
            if k not in control_keys
        }

    config["_resolved_init_payload"] = init_payload
    config["_worker_key"] = yaml.safe_dump(
        {
            "script_file": config["script_file"],
            "server_cwd": config.get("server_cwd"),
            "init_command": config["init_command"],
            "dock_command": config["dock_command"],
            "request_timeout": config.get("request_timeout"),
            "init_payload": init_payload,
        },
        sort_keys=True,
    )
    return config


def _prepare_ligand_payload(mol, config, ring_sample=False):
    payload_type = config["ligand_payload_type"]

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

    if payload_type == "mol_block":
        return [Chem.MolToMolBlock(mol)]

    raise ValueError(f"Unsupported ligand_payload_type: {payload_type}")


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
        command=config["script_file"],
        startup_timeout=config.get("startup_timeout", 120),
        request_timeout=config.get("request_timeout"),
        cwd=config.get("server_cwd"),
    )

    response = client.request(
        command=config["init_command"],
        payload=config["_resolved_init_payload"],
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


def _extract_batch_results(response, config):
    if not isinstance(response, dict):
        return []

    result_items_key = config["result_items_key"]
    payload = _response_payload(response)

    if isinstance(payload.get(result_items_key), list):
        return payload[result_items_key]

    logger.warning(
        f"Received response for 'payload' -> {result_items_key!r} is not a list. This is a break of API.")

    #
    # if isinstance(response.get(result_items_key), list):
    #     return response[result_items_key]
    #
    # score_key = config["score_key"]
    # pose_key = config["pose_key"]
    # mol_block_key = config["mol_block_key"]
    #
    # if score_key in payload and (pose_key in payload or mol_block_key in payload):
    #     return [payload]
    #
    # if score_key in response and (pose_key in response or mol_block_key in response):
    #     return [response]

    return []


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


def mol_dock(mol, config, ring_sample=False):
    mol_id = mol.GetProp("_Name")

    try:
        config = _parse_config(config)
        ligand_payload = _prepare_ligand_payload(mol, config, ring_sample=ring_sample)
        if ligand_payload is None:
            return mol_id, None
    except Exception as e:
        logger.warning("Ligand preparation failed for %s: %s", mol_id, e)
        return mol_id, None

    start_time = timeit.default_timer()

    try:
        client = _get_worker_client(config)
        response = client.request(
            command=config["dock_command"],
            payload={
                "molecule_id": mol_id,
                # "ring_sample": bool(ring_sample),
                config["ligand_payload_type"]: ligand_payload,
            },
            timeout=config.get("request_timeout"),
        )
        _check_response_ok(response, context="dock")
    except Exception as e:
        logger.warning("Docking request failed for %s: %s", mol_id, e)
        return mol_id, None

    results = _extract_batch_results(response, config)
    dock_output_conformer_list = []

    score_key = config["score_key"]
    pose_key = config["pose_key"]
    mol_block_key = config["mol_block_key"]

    raw_format = config.get("raw_format", "pdbqt")

    for item in results:
        docking_score = _safe_float(item.get(score_key))
        if docking_score is None:
            continue

        mol_block = item.get(mol_block_key)
        raw_block = item.get(pose_key)

        if isinstance(mol_block, str) and mol_block.strip():
            direct_output = {
                "docking_score": docking_score,
                "mol_block": mol_block,
            }
            if isinstance(raw_block, str):
                direct_output["raw_block"] = raw_block
            dock_output_conformer_list.append(direct_output)
            continue

        if not isinstance(raw_block, str) or not raw_block.strip():
            continue

        # parse raw_block if mol_block is absent
        if raw_format == "pdbqt":

            if "MODEL" not in raw_block:
                continue

            try:
                _, _pdbqt2molblock = _ensure_preparation_functions()
                mol_block = _pdbqt2molblock(raw_block.split("MODEL", 1)[1], mol, mol_id)
            except Exception:
                logger.exception("Failed to convert pdbqt pose for %s", mol_id)
                continue

            dock_output_conformer_list.append(
                {
                    "docking_score": docking_score,
                    "mol_block": mol_block,
                    "raw_block": raw_block,
                }
            )

        elif raw_format == "sdf":

            dock_output_conformer_list.append(
                {
                    "docking_score": docking_score,
                    "mol_block": raw_block.split('$$$$\n')[0],
                    "raw_block": raw_block,
                }
            )

        else:
            raise NotImplementedError(f'Output raw docking format {raw_format!r} is not implemented.')


    dock_time = round(timeit.default_timer() - start_time, 1)

    if dock_output_conformer_list:
        output = _choose_best(dock_output_conformer_list, config.get("score_mode", "min"))
        mol_name = _extract_molecule_name(response)
        if mol_name is not None and "mol_block" in output:
            output["mol_block"] = _set_mol_block_name(output["mol_block"], mol_name)
        output["dock_time"] = dock_time
    else:
        output = None

    return mol_id, output


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
