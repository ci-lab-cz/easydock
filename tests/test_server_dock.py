import sys
import textwrap

import yaml
from rdkit import Chem

import easydock.server_dock as server_dock


def _make_fake_server_script(tmp_path):
    script_path = tmp_path / "fake_docking_server.py"
    script_path.write_text(
        textwrap.dedent(
            """\
            import json
            import sys

            init_calls = 0
            dock_calls = 0

            for line in sys.stdin:
                if not line.strip():
                    continue
                req = json.loads(line)
                req_id = req.get("id")
                cmd = req.get("command")

                if cmd == "init":
                    init_calls += 1
                    resp = {"id": req_id, "status": "ok"}
                elif cmd == "dock_batch":
                    dock_calls += 1
                    ligands = req.get("payload", {}).get("ligands_pdbqt", [])
                    results = []
                    for i, _ in enumerate(ligands, 1):
                        results.append({
                            "docking_score": -1.0 * i,
                            "pdb_block": "MODEL 1\\nENDMDL\\n"
                        })
                    resp = {"id": req_id, "status": "ok", "results": results}
                elif cmd == "stats":
                    resp = {
                        "id": req_id,
                        "status": "ok",
                        "payload": {"init_calls": init_calls, "dock_calls": dock_calls},
                    }
                else:
                    resp = {"id": req_id, "status": "error", "error": "unknown command"}

                print(json.dumps(resp), flush=True)
            """
        ),
        encoding="utf-8",
    )
    return script_path


def _make_fake_smiles_server_script(tmp_path):
    script_path = tmp_path / "fake_smiles_server.py"
    script_path.write_text(
        textwrap.dedent(
            """\
            import json
            import sys

            init_calls = 0
            dock_calls = 0
            last_payload = {}

            for line in sys.stdin:
                if not line.strip():
                    continue
                req = json.loads(line)
                req_id = req.get("id")
                cmd = req.get("command")

                if cmd == "init":
                    init_calls += 1
                    resp = {"id": req_id, "status": "ok"}
                elif cmd == "dock_batch":
                    dock_calls += 1
                    last_payload = req.get("payload", {})
                    ligands = last_payload.get("ligands_smiles", [])
                    results = []
                    for _ in ligands:
                        results.append({
                            "docking_score": -7.3,
                            "mol_block": "mock_mol_block"
                        })
                    resp = {"id": req_id, "status": "ok", "results": results}
                elif cmd == "stats":
                    resp = {
                        "id": req_id,
                        "status": "ok",
                        "payload": {
                            "init_calls": init_calls,
                            "dock_calls": dock_calls,
                            "last_payload": last_payload,
                        },
                    }
                else:
                    resp = {"id": req_id, "status": "error", "error": "unknown command"}

                print(json.dumps(resp), flush=True)
            """
        ),
        encoding="utf-8",
    )
    return script_path


def test_server_backend_reuses_worker_client(tmp_path, monkeypatch):
    script_path = _make_fake_server_script(tmp_path)
    config_path = tmp_path / "server_config.yml"
    config = {
        "script_file": "{} {}".format(sys.executable, script_path),
        "startup_timeout": 5,
        "request_timeout": 5,
    }
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    monkeypatch.setattr(
        server_dock,
        "ligand_preparation",
        lambda mol, boron_replacement=True, ring_sample=False: ["MODEL 1\nENDMDL\n"],
    )
    monkeypatch.setattr(server_dock, "pdbqt2molblock", lambda pose, mol, mol_id: "mock_mol_block")

    mol1 = Chem.MolFromSmiles("CCO")
    mol1.SetProp("_Name", "mol_1_0")
    mol2 = Chem.MolFromSmiles("CCN")
    mol2.SetProp("_Name", "mol_2_0")

    server_dock._reset_worker_state_for_tests()
    try:
        mol_id1, out1 = server_dock.mol_dock(mol1, str(config_path))
        mol_id2, out2 = server_dock.mol_dock(mol2, str(config_path))

        assert mol_id1 == "mol_1_0"
        assert mol_id2 == "mol_2_0"
        assert out1 is not None
        assert out2 is not None
        assert out1["docking_score"] == -1.0
        assert out2["docking_score"] == -1.0

        parsed_cfg = server_dock._parse_config(str(config_path))
        client = server_dock._get_worker_client(parsed_cfg)
        stats_response = client.request("stats", payload={})
        stats_payload = stats_response.get("payload", stats_response)

        assert stats_payload["init_calls"] == 1
        assert stats_payload["dock_calls"] == 2
    finally:
        server_dock._reset_worker_state_for_tests()


def test_server_backend_supports_smiles_payload_and_direct_mol_block(tmp_path, monkeypatch):
    script_path = _make_fake_smiles_server_script(tmp_path)
    config_path = tmp_path / "smiles_server_config.yml"
    config = {
        "script_file": "{} {}".format(sys.executable, script_path),
        "startup_timeout": 5,
        "request_timeout": 5,
        "ligand_payload_type": "smiles",
        "ligand_payload_key": "ligands_smiles",
        "mol_block_key": "mol_block",
    }
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    def _unexpected_ligand_preparation(*_args, **_kwargs):
        raise AssertionError("ligand_preparation should not be called for smiles payload")

    monkeypatch.setattr(server_dock, "ligand_preparation", _unexpected_ligand_preparation)

    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", "mol_smiles_0")

    server_dock._reset_worker_state_for_tests()
    try:
        mol_id, output = server_dock.mol_dock(mol, str(config_path), ring_sample=True)
        assert mol_id == "mol_smiles_0"
        assert output is not None
        assert output["docking_score"] == -7.3
        assert output["mol_block"] == "mock_mol_block"
        assert "pdb_block" not in output

        parsed_cfg = server_dock._parse_config(str(config_path))
        client = server_dock._get_worker_client(parsed_cfg)
        stats_response = client.request("stats", payload={})
        stats_payload = stats_response.get("payload", stats_response)

        assert stats_payload["init_calls"] == 1
        assert stats_payload["dock_calls"] == 1
        assert stats_payload["last_payload"]["ligands_smiles"] == ["CCO"]
        assert stats_payload["last_payload"]["molecule_id"] == "mol_smiles_0"
        assert stats_payload["last_payload"]["ring_sample"] is True
    finally:
        server_dock._reset_worker_state_for_tests()
