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

