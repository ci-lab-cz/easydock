"""
Tests of the rDock server (containers/rdock) with a stub rbdock, and of the mol3d ligand format of server docking.
"""
import importlib.util
import logging
import os
import sqlite3
import stat
import sys
import textwrap

import pytest
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem

from easydock.dock.server_dock import _prepare_ligand_payload

RDOCK_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'containers', 'rdock')

STUB_RBDOCK = textwrap.dedent('''\
    #!/usr/bin/env python3
    """Stub rbdock: every input record is written 3 times to <root>.sd with SCORE/SCORE.INTER fields.
    SCORE decreases with the copy index, SCORE.INTER is lowest for the copy with the second-best SCORE."""
    import os, sys
    args = sys.argv[1:]
    opts = dict(zip(args[::2], args[1::2]))
    if 'FAIL' in os.getcwd():   # the server runs rbdock inside the receptor folder with a bare -r name
        print('RBT_FILE_READ_ERROR at src/exe/rbdock.cxx, line 368')
        print('Cavity file not found')
        sys.exit(0)
    if 'REJECT' in os.getcwd():  # a ligand without the features of the mandatory pharmacophore restraints
        print('RBT_LIGAND_ERROR at src/lib/RbtConstraint.cxx, line 297')
        print('The ligand has only 0 aromatic ring(s) (1 required)')
        print('END OF RUN')
        sys.exit(0)
    records = [r + '$$$$\\n' for r in open(opts['-i']).read().split('$$$$\\n') if r.strip()]
    with open(opts['-o'] + '.sd', 'w') as f:
        for rec in records:
            for i in range(3):
                score = 10.0 - i                # best SCORE (8.0) is copy 2
                inter = [5.0, 1.0, 3.0][i]      # lowest SCORE.INTER (1.0) is copy 1
                f.write(rec.replace('$$$$\\n', '') + f'>  <SCORE>\\n{score}\\n\\n>  <SCORE.INTER>\\n{inter}\\n\\n$$$$\\n')
    print('END OF RUN')
    ''')


@pytest.fixture(scope='module')
def server_module():
    sys.path.insert(0, RDOCK_DIR)
    spec = importlib.util.spec_from_file_location('rdock_server', os.path.join(RDOCK_DIR, 'rdock_server.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope='module')
def prep_module():
    spec = importlib.util.spec_from_file_location('rdock_prepare_receptor',
                                                  os.path.join(RDOCK_DIR, 'rdock_prepare_receptor.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def mol_block_3d(smi, name):
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):  # move away from the origin to check centring
        x, y, z = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (x + 50, y - 120, z + 7))
    mol.SetProp('_Name', name)
    return Chem.MolToMolBlock(mol)


def test_ligand_record_charges_title_and_centring(server_module):
    block = mol_block_3d('C[NH3+].[O-]C(=O)C', 'orig')
    record = server_module.ligand_record(block, 'lig_1')
    lines = record.split('\n')
    assert lines[0] == 'lig_1'
    assert record.endswith('M  END\n$$$$\n')
    n_atoms = int(lines[3][0:3])
    atom_lines = lines[4:4 + n_atoms]
    charge_codes = [int(ln[36:39]) for ln in atom_lines]
    symbols = [ln[31:34].strip() for ln in atom_lines]
    assert charge_codes.count(3) == 1 and symbols[charge_codes.index(3)] == 'N'   # +1
    assert charge_codes.count(5) == 1 and symbols[charge_codes.index(5)] == 'O'   # -1
    coords = [[float(ln[0:10]), float(ln[10:20]), float(ln[20:30])] for ln in atom_lines]
    for k in range(3):
        assert abs(sum(c[k] for c in coords) / n_atoms) < 1e-3
    assert Chem.MolFromMolBlock(record.replace('$$$$\n', ''), removeHs=False) is not None


def test_ligand_record_rejects_2d(server_module):
    mol = Chem.MolFromSmiles('CCO')
    with pytest.raises(ValueError, match='3D'):
        server_module.ligand_record(Chem.MolToMolBlock(mol), 'x')


def test_mol3d_payload():
    mol = Chem.MolFromSmiles('C1CCCCC1C(=O)O')
    mol.SetProp('_Name', 'acid')
    blocks = _prepare_ligand_payload(mol, {'ligand_in_format': 'mol3d'})
    assert len(blocks) == 1
    m = Chem.MolFromMolBlock(blocks[0], removeHs=False)
    assert m.GetConformer().Is3D()
    assert sum(a.GetAtomicNum() == 1 for a in m.GetAtoms()) == 12
    blocks = _prepare_ligand_payload(mol, {'ligand_in_format': 'mol3d'}, ring_sample=True)
    assert len(blocks) >= 1


@pytest.fixture
def stub_receptor(tmp_path, monkeypatch):
    stub = tmp_path / 'rbdock'
    stub.write_text(STUB_RBDOCK)
    stub.chmod(stub.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv('PATH', str(tmp_path) + os.pathsep + os.environ.get('PATH', ''))
    monkeypatch.delenv('RBT_ROOT', raising=False)
    rec = tmp_path / 'rec'
    rec.mkdir()
    (rec / 'receptor.prm').write_text('RBT_PARAMETER_FILE_V1.00\nRECEPTOR_FILE protein.mol2\n')
    (rec / 'receptor.as').write_bytes(b'cavity')
    return rec


def test_dock_with_stub_rbdock(server_module, stub_receptor):
    server = server_module._RDockServer()
    res = server.init({'receptor': str(stub_receptor / 'receptor.prm'), 'n_poses': 2, 'seed': None})
    assert res == {'status': 'ok'}, res
    payload = {'lig_a': mol_block_3d('CCO', 'a'), 'lig_b': mol_block_3d('c1ccccc1C(=O)[O-]', 'b'), 'bad': 'not a mol'}
    res = server.dock(payload)
    assert res['status'] == 'ok'
    results = res['results']
    assert results['bad'] is None
    for name in ('lig_a', 'lig_b'):
        out = results[name]
        assert out['docking_score'] == 3.0           # SCORE.INTER of the pose with the best SCORE (copy 2)
        assert out['raw_block'].count('$$$$\n') == 2  # n_poses
        assert server_module.sd_field(out['raw_block'].split('$$$$\n')[0], 'SCORE') == '8.0'
        assert out['mol_block'].split('\n')[0] == name
        mol = Chem.MolFromMolBlock(out['mol_block'], removeHs=False)
        assert mol is not None and mol.GetConformer().Is3D()
    assert Chem.MolToSmiles(Chem.RemoveHs(Chem.MolFromMolBlock(results['lig_b']['mol_block']))) == \
        Chem.MolToSmiles(Chem.MolFromSmiles('c1ccccc1C(=O)[O-]'))

    # one field for both ranking and score
    server.init({'receptor': str(stub_receptor / 'receptor.prm'), 'sort_field': 'SCORE.INTER',
                 'score_field': 'SCORE.INTER'})
    out = server.dock({'lig_a': payload['lig_a']})['results']['lig_a']
    assert out['docking_score'] == 1.0


def test_dock_fatal_error_fails_request(server_module, stub_receptor, tmp_path):
    fail = tmp_path / 'FAIL'
    fail.mkdir()
    (fail / 'receptor.prm').write_text('RBT_PARAMETER_FILE_V1.00\n')
    (fail / 'receptor.as').write_bytes(b'cavity')
    server = server_module._RDockServer()
    assert server.init({'receptor': str(fail / 'receptor.prm')})['status'] == 'ok'
    res = server.dock({'lig_a': mol_block_3d('CCO', 'a')})
    assert res['status'] == 'error' and 'RBT_FILE_READ_ERROR' in res['error']


def test_init_errors(server_module, stub_receptor, tmp_path):
    server = server_module._RDockServer()
    assert 'required' in server.init({})['error']
    (tmp_path / 'nocav.prm').write_text('RBT_PARAMETER_FILE_V1.00\n')
    assert 'cavity' in server.init({'receptor': str(tmp_path / 'nocav.prm')})['error']
    assert server.dock({'x': 'y'})['status'] == 'error'


PDB = textwrap.dedent('''\
    ATOM      1  CG  HIS A   5      10.000  10.000  10.000  1.00  0.00           C
    ATOM      2  ND1 HIS A   5      11.000  10.000  10.000  1.00  0.00           N
    ATOM      3  HD1 HIS A   5      11.500  10.500  10.000  1.00  0.00           H
    ATOM      4  N   GLY B   5      20.000  10.000  10.000  1.00  0.00           N
    HETATM    5  O   HOH A 100      30.000  10.000  10.000  1.00  0.00           O
    ''')

RAW_MOL2 = textwrap.dedent('''\
    @<TRIPOS>MOLECULE
    clean.pdb
     4 2 0 0 0
    SMALL
    GASTEIGER

    @<TRIPOS>ATOM
          1  CG         10.0000   10.0000   10.0000 C.3       5  HIS5        0.1000
          2  ND1        11.0000   10.0000   10.0000 N.pl3     5  HIS5       -0.2000
          3  HD1        11.5000   10.5000   10.0000 H         5  HIS5        0.3000
          4  N          20.0000   10.0000   10.0000 N.3       5  GLY5       -0.4000
    @<TRIPOS>BOND
         1     2     1    1
         2     3     2    1
    ''')


def test_mol2_fixup(prep_module, tmp_path):
    pdb = tmp_path / 'p.pdb'
    pdb.write_text(PDB)
    atoms = prep_module.read_pdb_atoms(str(pdb))
    assert len(atoms) == 4 and prep_module.has_hydrogens(atoms)      # the water is dropped
    raw = tmp_path / 'raw.mol2'
    raw.write_text(RAW_MOL2)
    out = tmp_path / 'protein.mol2'
    assert prep_module.mol2_fixup(atoms, str(raw), str(out)) == 2
    lines = out.read_text().split('\n')
    assert lines[2] == ' 4 2 2 0 0'
    atom_lines = lines[lines.index('@<TRIPOS>ATOM') + 1:lines.index('@<TRIPOS>BOND')]
    fields = [ln.split() for ln in atom_lines]
    assert [f[5] for f in fields] == ['C.ar', 'N.ar', 'H', 'N.3']   # ring atoms typed by name
    assert [f[6] for f in fields] == ['1', '1', '1', '2']           # one subst_id per residue, chains apart
    assert [f[7] for f in fields] == ['HIS5', 'HIS5', 'HIS5', 'GLY5']
    assert [f[8] for f in fields] == ['0.1000', '-0.2000', '0.3000', '-0.4000']
    bonds = lines[lines.index('@<TRIPOS>BOND') + 1:lines.index('@<TRIPOS>SUBSTRUCTURE')]
    assert [b.split() for b in bonds] == [['1', '1', '2', '1'], ['2', '2', '3', '1']]  # low atom id first
    subs = [s.split() for s in lines[lines.index('@<TRIPOS>SUBSTRUCTURE') + 1:] if s.strip()]
    assert [(s[0], s[1], s[2], s[5], s[6]) for s in subs] == [('1', 'HIS5', '1', 'A', 'HIS'),
                                                              ('2', 'GLY5', '4', 'B', 'GLY')]


def test_reference_ligand_is_kekulized(prep_module, tmp_path):
    mol = Chem.AddHs(Chem.MolFromSmiles('c1ccccc1O'))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    block = Chem.MolToMolBlock(mol, kekulize=False)   # aromatic bonds as type 4, as Marvin writes them
    assert any(ln.endswith('  4  0') or ln[6:9] == '  4' for ln in block.split('\n')[4 + mol.GetNumAtoms():])
    src = tmp_path / 'ref.sdf'
    src.write_text(block + '$$$$\n')
    dst = tmp_path / 'ref_ligand.sdf'
    prep_module.write_reference_ligand(str(src), str(dst))
    out = Chem.MolFromMolFile(str(dst), removeHs=False)
    assert out.GetNumAtoms() == mol.GetNumAtoms() and out.GetConformer().Is3D()
    bond_lines = dst.read_text().split('\n')[4 + mol.GetNumAtoms():4 + mol.GetNumAtoms() + mol.GetNumBonds()]
    assert all(int(ln[6:9]) in (1, 2, 3) for ln in bond_lines)


def test_bad_log_pattern(prep_module):
    log = 'VOL_INCR\t0\n\nRBT_MODEL_ERROR at src/lib/RbtMdlFileSource.cxx, line 280\nC32 makes a bond\n'
    assert prep_module.BAD_LOG.search(log).group(0) == \
        'RBT_MODEL_ERROR at src/lib/RbtMdlFileSource.cxx, line 280\nC32 makes a bond'
    assert prep_module.BAD_LOG.search('**WARNING** RBT_MODEL_ERROR at src/x.cxx, line 393\nN16 makes too few bonds\n') \
        is None


def test_mapper_and_prm(prep_module, tmp_path):
    lines = prep_module.mapper_section(center=(1, 2.5, -3), radius=12)
    assert '    CENTER (1.000,2.500,-3.000)' in lines and '    RADIUS 12.0' in lines
    assert '    LARGE_SPHERE 4.0' in lines
    lines = prep_module.mapper_section(reflig_name='ref_ligand.sdf')
    assert '    REF_MOL ref_ligand.sdf' in lines and '    RADIUS 6.0' in lines
    assert not any('LARGE_SPHERE' in ln for ln in lines)  # rejected by the ligand site mapper
    prm = tmp_path / 'r.prm'
    prep_module.write_prm(str(prm), lines, receptor_flex=3.0)
    text = prm.read_text().split('\n')
    assert text[0] == 'RBT_PARAMETER_FILE_V1.00'
    assert 'RECEPTOR_FLEX 3.0' in text and text.index('RECEPTOR_FLEX 3.0') < text.index('SECTION MAPPER')
    with pytest.raises(ValueError):
        prep_module.mapper_section()


def test_dock_ligand_rejection_is_not_fatal(server_module, stub_receptor, tmp_path, caplog):
    rej = tmp_path / 'REJECT'
    rej.mkdir()
    (rej / 'receptor.prm').write_text('RBT_PARAMETER_FILE_V1.00\n')
    (rej / 'receptor.as').write_bytes(b'cavity')
    server = server_module._RDockServer()
    assert server.init({'receptor': str(rej / 'receptor.prm')})['status'] == 'ok'
    with caplog.at_level(logging.WARNING):
        res = server.dock({'lig_a': mol_block_3d('CCO', 'a')})
    assert res == {'status': 'ok', 'results': {'lig_a': None}}
    assert 'rejected lig_a' in caplog.text and 'aromatic ring' in caplog.text


def test_parse_restraints(prep_module, tmp_path):
    f = tmp_path / 'ph4.const'
    f.write_text('# a comment line\n1.0, 2.0, 3.0, 1.5, Acc\n\n-1 0 0.5 2 Don  # the donor hydrogen\n')
    assert prep_module.parse_restraints(str(f)) == [(1.0, 2.0, 3.0, 1.5, 'Acc'), (-1.0, 0.0, 0.5, 2.0, 'Don')]
    assert prep_module.parse_restraints([[1, 2, 3, 1.5, 'Aro'], '4 5 6 1 Hyd']) == \
        [(1.0, 2.0, 3.0, 1.5, 'Aro'), (4.0, 5.0, 6.0, 1.0, 'Hyd')]
    assert prep_module.parse_restraints('1 2 3 1 Any\n4 5 6 1 Cat\n') == [(1.0, 2.0, 3.0, 1.0, 'Any'),
                                                                        (4.0, 5.0, 6.0, 1.0, 'Cat')]
    for bad in ('missing.const', [[1, 2, 3, 1.5, 'Acceptor']], [[1, 2, 3, 1.5, 'acc']], [[1, 2, 3, 'Acc']],
                [[1, 2, 3, 0, 'Acc']], [], '# nothing\n\n', 'x y z 1 Acc\n', 5):
        with pytest.raises(ValueError):
            prep_module.parse_restraints(bad)


def test_pharma_section_and_prm(prep_module, tmp_path):
    out = tmp_path / 'rec'
    out.mkdir()
    lines = prep_module.prepare_pharma({'constraints_file': [[1, 2, 3, 1.5, 'Acc']],
                                        'optional_file': '0 0 0 1 Don\n1 1 1 1 Aro\n', 'nopt': 1, 'weight': 2},
                                       str(out))
    assert (out / 'pharma_mandatory.const').read_text() == '1.000 2.000 3.000 1.50 Acc\n'
    assert (out / 'pharma_optional.const').read_text().splitlines() == ['0.000 0.000 0.000 1.00 Don',
                                                                        '1.000 1.000 1.000 1.00 Aro']
    assert lines == ['SECTION PHARMA', '    SCORING_FUNCTION RbtPharmaSF', '    WEIGHT 2.0',
                     '    CONSTRAINTS_FILE pharma_mandatory.const', '    OPTIONAL_FILE pharma_optional.const',
                     '    NOPT 1', 'END_SECTION']
    prm = tmp_path / 'r.prm'
    prep_module.write_prm(str(prm), prep_module.mapper_section(center=(0, 0, 0)), pharma_lines=lines)
    text = prm.read_text().split('\n')
    assert text.index('SECTION PHARMA') > text.index('SECTION CAVITY') and 'END_SECTION' == text[-2]

    lines = prep_module.prepare_pharma({'constraints_file': '1 2 3 1 Acc\n', 'write_errors': True}, str(out))
    assert lines == ['SECTION PHARMA', '    SCORING_FUNCTION RbtPharmaSF', '    WEIGHT 1.0',
                     '    CONSTRAINTS_FILE pharma_mandatory.const', '    WRITE_ERRORS TRUE', 'END_SECTION']
    for bad in ({'optional_file': '1 2 3 1 Acc\n'},                                   # mandatory file missing
                {'constraints_file': '1 2 3 1 Acc\n', 'nopt': 1},                      # nopt without optional points
                {'constraints_file': '1 2 3 1 Acc\n', 'optional_file': '1 2 3 1 Acc\n', 'nopt': 2},
                {'constraints_file': '1 2 3 1 Acc\n', 'weight': -1},
                {'constraints_file': '1 2 3 1 Acc\n', 'constraints': 'x'},             # unknown key
                'pharma.const'):
        with pytest.raises(ValueError):
            prep_module.prepare_pharma(bad, str(out))


def test_init_pharma(server_module, stub_receptor, tmp_path, monkeypatch):
    server = server_module._RDockServer()
    pharma = {'constraints_file': [[0, 0, 0, 1, 'Any']]}
    res = server.init({'receptor': str(stub_receptor / 'receptor.prm'), 'pharma': pharma})
    assert res['status'] == 'error' and 'SECTION PHARMA' in res['error']
    assert 'mapping' in server.init({'protein': 'p.pdb', 'pharma': 'pharma.const'})['error']

    calls = {}

    def fake_prepare(protein, out_dir, **kwargs):
        calls.update(kwargs)
        with open(os.path.join(out_dir, 'receptor.prm'), 'w') as f:
            f.write('RBT_PARAMETER_FILE_V1.00\n')
        with open(os.path.join(out_dir, 'receptor.as'), 'wb') as f:
            f.write(b'cavity')
        return {'prm': os.path.join(out_dir, 'receptor.prm'), 'cavity_volume': 1.0, 'n_residues': 1}

    monkeypatch.setattr(server_module.receptor_prep, 'prepare', fake_prepare)
    res = server.init({'protein': str(tmp_path / 'p.pdb'), 'center': [0, 0, 0], 'pharma': pharma,
                       'tmpdir': str(tmp_path)})
    assert res == {'status': 'ok'}, res
    assert calls['pharma'] == pharma and calls['center'] == [0, 0, 0]
    assert server.dock({'lig_a': mol_block_3d('CCO', 'a')})['results']['lig_a']['docking_score'] == 3.0
    server.cleanup()


def test_pharma_files_restored_on_resume(tmp_path):
    """Restraint files given as paths are stored in the setup table like other config files and come back as temp
    files on resume; inline lists stay in the config text."""
    from easydock.session import (create_setup_table, save_session_args, save_session_config, restore_session,
                                  DB_USER_VERSION)
    const = tmp_path / 'mandatory.const'
    const.write_text('1 2 3 1.5 Acc\n')
    cfg = tmp_path / 'config.yml'
    cfg.write_text(yaml.safe_dump({'script_file': 'rdock.sif',
                                   'init_server': {'protein': 'p.pdb',
                                                   'pharma': {'constraints_file': str(const),
                                                              'optional_file': [[0, 0, 0, 1, 'Don']],
                                                              'nopt': 1}}}))
    db = tmp_path / 'session.db'
    with sqlite3.connect(db) as conn:
        create_setup_table(conn)
        save_session_args(conn, {'config': str(cfg)})
        save_session_config(conn, str(cfg))
        conn.execute(f'PRAGMA user_version = {DB_USER_VERSION}')
        keys = sorted(r[0] for r in conn.execute('SELECT key FROM setup'))
    assert keys == ['args', 'config', 'file:init_server.pharma.constraints_file']
    const.unlink()  # the original file is not needed on resume
    args, tmpfiles = restore_session(str(db), tmpdir=str(tmp_path))
    try:
        with open(args['config']) as f:
            restored = yaml.safe_load(f)['init_server']['pharma']
        assert restored['constraints_file'] != str(const) and restored['constraints_file'].endswith('.const')
        with open(restored['constraints_file']) as f:
            assert f.read() == '1 2 3 1.5 Acc\n'
        assert restored['optional_file'] == [[0, 0, 0, 1, 'Don']] and restored['nopt'] == 1
    finally:
        for f in tmpfiles:
            os.unlink(f)
