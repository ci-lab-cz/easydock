"""
Dummy docking program for generic docking tests.
Reads SMILES from stdin and writes 3 embedded conformers to stdout as poses.
Output format (sdf or pdb) is taken from the DUMMY_OUT_FORMAT environment variable.
"""
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

N_POSES = 3

out_format = os.environ.get('DUMMY_OUT_FORMAT', 'sdf')
mol = Chem.AddHs(Chem.MolFromSmiles(sys.stdin.read().strip()))
cids = AllChem.EmbedMultipleConfs(mol, N_POSES, randomSeed=42)

output = ''
for i, cid in enumerate(cids, 1):
    score = -10 + i
    if out_format == 'sdf':
        output += Chem.MolToMolBlock(mol, confId=cid) + f'>  <score>\n{score}\n\n$$$$\n'
    else:
        pdb_block = Chem.MolToPDBBlock(Chem.RemoveHs(mol), confId=cid).replace('END\n', '')
        output += f'MODEL {i}\nREMARK SCORE {score}\n{pdb_block}ENDMDL\n'
sys.stdout.write(output)
