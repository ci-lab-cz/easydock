"""
Tests of storing raw_format in a database and extracting poses from raw_block
for generic and server docking (issue #74).
"""
import os
import sqlite3
import subprocess
import sys
from argparse import Namespace

import pytest
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem

from easydock.database import (_get_raw_format, create_db, get_poses_from_raw_block, get_variables,
                               populate_setup_db)
from easydock.dock.generic_dock import _build_mol_block

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
DUMMY_GENERIC_DOCK = os.path.join(DATA_DIR, 'dummy_generic_dock.py')
DUMMY_SERVER = os.path.join(DATA_DIR, 'dummy_server.py')

INPUT_SMI = 'CC(=O)Nc1ccc(O)cc1\tparacetamol\nO=C(O)c1ccccc1OC(C)=O\taspirin\n'
EXPECTED_NAMES = [f'{name}:{i}' for name in ('paracetamol', 'aspirin') for i in (1, 2, 3)]


def write_yaml(path, data):
    with open(path, 'w') as f:
        yaml.safe_dump(data, f)
    return str(path)


def get_raw_format_from_db(db_fname):
    with sqlite3.connect(db_fname) as conn:
        return get_variables(conn, 'database', ['raw_format'])['raw_format']


@pytest.fixture
def mol_with_confs():
    mol = Chem.AddHs(Chem.MolFromSmiles('O=C(O)c1ccccc1N'))
    cids = list(AllChem.EmbedMultipleConfs(mol, 3, randomSeed=1))
    return mol, cids


# ---------------------------------------------------------------------------
# _get_raw_format
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('program, config, expected', [
    ('generic', {'ligand_out_format': 'sdf'}, 'sdf'),
    ('generic', {'ligand_out_format': 'pdb'}, 'pdb'),
    ('generic', {'ligand_out_format': 'pdbqt'}, 'pdbqt'),
    ('server', {'ligand_out_format': 'sdf'}, 'sdf'),
    ('vina', {'ligand_out_format': 'sdf'}, None),
    ('vina', {'raw_format': 'sdf'}, None),
    (None, {'ligand_out_format': 'sdf'}, None),
])
def test_get_raw_format(tmp_path, program, config, expected):
    assert _get_raw_format(write_yaml(tmp_path / 'config.yml', config), program) == expected


def test_get_raw_format_server_queries_info(tmp_path, monkeypatch):
    monkeypatch.setenv('DUMMY_SERVER_LOG', str(tmp_path / 'server.log'))
    config = write_yaml(tmp_path / 'config.yml', {'script_file': f'{sys.executable} {DUMMY_SERVER}'})
    assert _get_raw_format(config, 'server') == 'sdf'


def test_get_raw_format_invalid(tmp_path):
    config = write_yaml(tmp_path / 'config.yml', {'ligand_out_format': 'mol2'})
    with pytest.raises(ValueError):
        _get_raw_format(config, 'generic')


# ---------------------------------------------------------------------------
# populate_setup_db
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('args_extra, config, expected', [
    ({'program': 'generic'}, {'ligand_out_format': 'sdf'}, 'sdf'),
    ({'program': 'vina'}, {'exhaustiveness': 8}, 'pdbqt'),
    ({'program': 'vina'}, {'raw_format': 'sdf'}, 'sdf'),
    ({}, {'ligand_out_format': 'sdf'}, 'pdbqt'),  # args without program, e.g. from external tools
])
def test_populate_setup_db(tmp_path, args_extra, config, expected):
    db_fname = str(tmp_path / 'test.db')
    args = Namespace(config=write_yaml(tmp_path / 'config.yml', config), **args_extra)
    create_db(db_fname, args)
    populate_setup_db(db_fname, args)
    assert get_raw_format_from_db(db_fname) == expected


# ---------------------------------------------------------------------------
# get_poses_from_raw_block
# ---------------------------------------------------------------------------

def test_get_poses_from_raw_block_sdf(mol_with_confs):
    mol, cids = mol_with_confs
    raw_block = ''.join(Chem.MolToMolBlock(mol, confId=cid) + '$$$$\n' for cid in cids)
    poses = get_poses_from_raw_block(raw_block, 'sdf', mol, 'mol_0', [2, 3])
    assert [i for i, _ in poses] == [2, 3]
    for i, block in poses:
        m = Chem.MolFromMolBlock(block)
        assert m.GetProp('_Name') == f'mol_0:{i}'
        assert Chem.MolToSmiles(m) == Chem.MolToSmiles(Chem.RemoveHs(mol))


def test_get_poses_from_raw_block_pdb(mol_with_confs):
    mol, cids = mol_with_confs
    mol_noh = Chem.RemoveHs(mol)
    raw_block = ''.join(f'MODEL {i}\n' + Chem.MolToPDBBlock(mol_noh, confId=cid).replace('END\n', '') + 'ENDMDL\n'
                        for i, cid in enumerate(cids, 1))
    template = Chem.MolFromMolBlock(_build_mol_block(raw_block, 'pdb', mol, 'mol_0'))
    assert template.GetProp('_Name') == 'mol_0'
    assert Chem.MolToSmiles(template) == Chem.MolToSmiles(mol_noh)

    poses = get_poses_from_raw_block(raw_block, 'pdb', template, 'mol_0', [2, 3, 4])  # pose 4 does not exist
    assert [i for i, _ in poses] == [2, 3]
    for i, block in poses:
        m = Chem.MolFromMolBlock(block)
        assert m.GetProp('_Name') == f'mol_0:{i}'
        assert Chem.MolToSmiles(m) == Chem.MolToSmiles(mol_noh)
        ref_pos = mol_noh.GetConformer(cids[i - 1]).GetPositions()
        assert abs(m.GetConformer().GetPositions() - ref_pos).max() < 0.01


def test_get_poses_from_raw_block_unknown_format(mol_with_confs):
    mol, _ = mol_with_confs
    assert get_poses_from_raw_block('anything', 'xyz', mol, 'mol_0', [2]) == []


def test_build_mol_block_sdf_title(mol_with_confs):
    mol, cids = mol_with_confs
    raw_block = Chem.MolToMolBlock(mol, confId=cids[0]) + '>  <score>\n-9\n\n$$$$\n'
    mol_block = _build_mol_block(raw_block, 'sdf', mol, 'mol_0')
    assert mol_block.split('\n', 1)[0] == 'mol_0'
    assert '<score>' not in mol_block
    assert Chem.MolFromMolBlock(mol_block) is not None


# ---------------------------------------------------------------------------
# end-to-end: run_dock + get_sdf_from_dock_db
# ---------------------------------------------------------------------------

def run_docking(tmp_path, program, config, env_extra):
    input_fname = tmp_path / 'input.smi'
    input_fname.write_text(INPUT_SMI)
    db_fname = str(tmp_path / 'output.db')
    sdf_fname = str(tmp_path / 'output.sdf')
    env = dict(os.environ, **env_extra)
    subprocess.run([sys.executable, '-m', 'easydock.run_dock', '-i', str(input_fname), '-o', db_fname,
                    '--program', program, '--config', config, '-c', '1'],
                   check=True, env=env, capture_output=True)
    subprocess.run([sys.executable, '-m', 'easydock.get_sdf_from_dock_db', '-i', db_fname, '-o', sdf_fname,
                    '--poses', '1', '2', '3'],
                   check=True, env=env, capture_output=True)
    names = [m.GetProp('_Name') for m in Chem.SDMolSupplier(sdf_fname)]
    return get_raw_format_from_db(db_fname), names


@pytest.mark.parametrize('out_format', ['sdf', 'pdb'])
def test_generic_docking_poses(tmp_path, out_format):
    config = write_yaml(tmp_path / 'config.yml', {
        'script_file': DUMMY_GENERIC_DOCK,
        'ligand_in_format': 'smiles',
        'ligand_out_format': out_format,
        'parse_score': {'score_field': 'score', 'score_regex': r'REMARK SCORE\s+(-?[\d.]+)'},
    })
    raw_format, names = run_docking(tmp_path, 'generic', config, {'DUMMY_OUT_FORMAT': out_format})
    assert raw_format == out_format
    assert names == EXPECTED_NAMES


@pytest.mark.parametrize('ligand_out_format, expected_starts', [
    (None, 3),   # info is queried in populate_setup_db and run_dock, plus a docking worker
    ('sdf', 2),  # info is queried only in run_dock, plus a docking worker
])
def test_server_docking_poses(tmp_path, ligand_out_format, expected_starts):
    config = {'script_file': f'{sys.executable} {DUMMY_SERVER}', 'init_server': {}}
    if ligand_out_format:
        config['ligand_out_format'] = ligand_out_format
    config = write_yaml(tmp_path / 'config.yml', config)
    log_fname = tmp_path / 'server.log'
    raw_format, names = run_docking(tmp_path, 'server', config, {'DUMMY_SERVER_LOG': str(log_fname)})
    assert raw_format == 'sdf'
    assert names == EXPECTED_NAMES
    assert len(log_fname.read_text().splitlines()) == expected_starts
