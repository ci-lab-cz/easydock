"""
Dummy docking server for server docking tests (JSON Lines over STDIN/STDOUT).
Returns 3 embedded conformers per molecule as SDF poses.
Each start of the server is logged to a file given by the DUMMY_SERVER_LOG environment variable.
"""
import json
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

N_POSES = 3

log_fname = os.environ.get('DUMMY_SERVER_LOG')
if log_fname:
    with open(log_fname, 'a') as f:
        f.write(f'start {os.getpid()}\n')


def dock(name, smi):
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    cids = AllChem.EmbedMultipleConfs(mol, N_POSES, randomSeed=42)
    raw_block = ''.join(f'{name}\n' + Chem.MolToMolBlock(mol, confId=cid).split('\n', 1)[1] + '$$$$\n'
                        for cid in cids)
    return {'docking_score': -9.0, 'raw_block': raw_block}


for line in sys.stdin:
    request = json.loads(line)
    command, payload = request['command'], request.get('payload', {})
    if command == 'info':
        response = {'ligand_in_format': 'smiles', 'ligand_out_format': 'sdf', 'score_mode': 'min', 'batch_size': 2}
    elif command == 'init':
        response = {'status': 'ok'}
    elif command == 'dock':
        response = {'status': 'ok', 'results': {name: dock(name, smi) for name, smi in payload.items()}}
    else:
        response = {'status': 'error', 'error': f'Unknown command: {command}'}
    sys.stdout.write(json.dumps({'id': request['id'], 'payload': response}) + '\n')
    sys.stdout.flush()
