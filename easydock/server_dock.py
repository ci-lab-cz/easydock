#!/usr/bin/env python3

import logging
import timeit
import traceback
from functools import lru_cache
from threading import Lock
from typing import Any, Dict, List, Optional

import yaml

from easydock.auxiliary import expand_path
from easydock.persistent_client import JsonLineProcessClient
from easydock.preparation_for_docking import ligand_preparation, pdbqt2molblock


_worker_lock = Lock()
_worker_client = None  # type: Optional[JsonLineProcessClient]
_worker_key = None  # type: Optional[str]

_control_config_keys = {
    "script_file",
    "init_payload",
    "init_command",
    "dock_command",
    "result_items_key",
    "score_key",
    "pose_key",
    "score_mode",
    "startup_timeout",
    "request_timeout",
    "boron_replacement",
    "server_cwd",
}


def _safe_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _normalize_path_values(config: Dict[str, Any]) -> Dict[str, Any]:
    path_keys = ("script_file", "protein", "protein_setup", "server_cwd")
    for key in path_keys:
        value = config.get(key)
        if isinstance(value, str):
            config[key] = expand_path(value)
    return config


@lru_cache(maxsize=64)
def _parse_config(config_fname: str) -> Dict[str, Any]:
    with open(config_fname) as f:
        config = yaml.safe_load(f) or {}

    if "script_file" not in config:
        raise KeyError("Configuration for --program server must contain `script_file`")

    config = _normalize_path_values(config)
    config.setdefault("init_command", "init")
    config.setdefault("dock_command", "dock_batch")
    config.setdefault("result_items_key", "results")
    config.setdefault("score_key", "docking_score")
    config.setdefault("pose_key", "pdb_block")
    config.setdefault("score_mode", "min")
    config.setdefault("startup_timeout", 120)
    config.setdefault("request_timeout", None)
    config.setdefault("boron_replacement", True)

    if config["request_timeout"] is not None:
        config["request_timeout"] = _safe_float(config["request_timeout"])
    config["startup_timeout"] = _safe_float(config["startup_timeout"]) or 120.0

    init_payload = config.get("init_payload")
    if not isinstance(init_payload, dict):
        init_payload = {
            key: value
            for key, value in config.items()
            if key not in _control_config_keys
        }
    config["_resolved_init_payload"] = init_payload

    # Key used to decide whether existing client can be reused in a worker process.
    config["_worker_key"] = yaml.safe_dump(
        {
            "script_file": config["script_file"],
            "init_payload": config["_resolved_init_payload"],
            "init_command": config["init_command"],
            "dock_command": config["dock_command"],
            "request_timeout": config["request_timeout"],
            "server_cwd": config.get("server_cwd"),
        },
        sort_keys=True,
    )
    return config


def _response_payload(response: Dict[str, Any]) -> Dict[str, Any]:
    payload = response.get("payload")
    if isinstance(payload, dict):
        return payload
    return response


def _check_response_ok(response: Dict[str, Any], context: str) -> None:
    payload = _response_payload(response)
    status = payload.get("status", response.get("status", "ok"))
    if status not in ("ok", "OK", True, None):
        error_message = payload.get("error", response.get("error", "<no error message>"))
        raise RuntimeError("{} failed: {}".format(context, error_message))


def _create_client(config: Dict[str, Any]) -> JsonLineProcessClient:
    client = JsonLineProcessClient(
        command=config["script_file"],
        startup_timeout=float(config["startup_timeout"]),
        request_timeout=config["request_timeout"],
        cwd=config.get("server_cwd"),
    )
    init_response = client.request(
        config["init_command"],
        payload=config["_resolved_init_payload"],
        timeout=config["request_timeout"],
    )
    _check_response_ok(init_response, "Docking server initialization")
    return client


def _get_worker_client(config: Dict[str, Any]) -> JsonLineProcessClient:
    global _worker_client, _worker_key
    with _worker_lock:
        must_create = (
            _worker_client is None
            or not _worker_client.is_alive()
            or _worker_key != config["_worker_key"]
        )
        if must_create:
            if _worker_client is not None:
                _worker_client.close()
            _worker_client = _create_client(config)
            _worker_key = config["_worker_key"]
        return _worker_client


def _extract_batch_results(response: Dict[str, Any], config: Dict[str, Any]) -> List[Dict[str, Any]]:
    payload = _response_payload(response)
    items_key = config["result_items_key"]
    if isinstance(payload.get(items_key), list):
        return payload[items_key]
    if isinstance(response.get(items_key), list):
        return response[items_key]

    score_key = config["score_key"]
    pose_key = config["pose_key"]
    if score_key in payload and pose_key in payload:
        return [payload]
    if score_key in response and pose_key in response:
        return [response]
    return []


def _choose_best(items: List[Dict[str, Any]], score_mode: str) -> Optional[Dict[str, Any]]:
    if not items:
        return None
    if score_mode == "max":
        return max(items, key=lambda x: x["docking_score"])
    return min(items, key=lambda x: x["docking_score"])


def mol_dock(mol, config, ring_sample=False):
    """
    Dock a molecule using a long-lived server process.

    The server is started once per worker process and reused for subsequent molecules.
    """
    config = _parse_config(config)
    mol_id = mol.GetProp("_Name")
    ligand_pdbqt_list = ligand_preparation(
        mol,
        boron_replacement=bool(config["boron_replacement"]),
        ring_sample=ring_sample,
    )

    if ligand_pdbqt_list is None:
        return mol_id, None

    start_time = timeit.default_timer()
    dock_output_conformer_list = []
    score_key = config["score_key"]
    pose_key = config["pose_key"]

    try:
        client = _get_worker_client(config)
        response = client.request(
            config["dock_command"],
            payload={
                "molecule_id": mol_id,
                "ring_sample": bool(ring_sample),
                "ligands_pdbqt": ligand_pdbqt_list,
            },
            timeout=config["request_timeout"],
        )
        _check_response_ok(response, "Docking request")
        response_items = _extract_batch_results(response, config)
    except Exception:
        logging.warning(
            "(server) Error caused by docking of %s\n%s",
            mol_id,
            traceback.format_exc(),
        )
        return mol_id, None

    for item in response_items:
        try:
            score = _safe_float(item.get(score_key))
            pdbqt_out = item.get(pose_key)
            if score is None or not isinstance(pdbqt_out, str):
                continue
            if "MODEL" not in pdbqt_out:
                logging.debug("(server) Missing MODEL section in docking output of %s", mol_id)
                continue
            mol_block = pdbqt2molblock(pdbqt_out.split("MODEL", 1)[1], mol, mol_id)
            dock_output_conformer_list.append(
                {
                    "docking_score": score,
                    "pdb_block": pdbqt_out,
                    "mol_block": mol_block,
                }
            )
        except Exception:
            logging.warning(
                "(server) Failed to parse one conformer for %s\n%s",
                mol_id,
                traceback.format_exc(),
            )

    dock_time = round(timeit.default_timer() - start_time, 1)
    logging.debug("(server) %s, docked nconf %s", mol_id, len(dock_output_conformer_list))

    output = _choose_best(dock_output_conformer_list, config["score_mode"])
    if output is None:
        return mol_id, None
    output["dock_time"] = dock_time
    return mol_id, output


def _reset_worker_state_for_tests() -> None:
    global _worker_client, _worker_key
    with _worker_lock:
        if _worker_client is not None:
            _worker_client.close()
        _worker_client = None
        _worker_key = None
