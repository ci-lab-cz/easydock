import subprocess
import os.path
import sys

from rdkit import Chem
from dimorphite_dl.dimorphite_dl import run as dimorphite_run

"""
Each protonation program should have two implemented functions:
1) protonate_xxx with two necessary arguments input_fname and output_fname (input and output filenames). 
Input file will be SMILES with names separated by tab. Output file format may be any (output file does not have 
an extension).
2) read_protonate_xxx, which takes a file name with protonated molecules in the format corresponding to 
the protonate_xxx function and returns a generator of tuples (SMILES, mol_name).
"""


def protonate_chemaxon(input_fname, output_fname, tautomerize=True):
    cmd_run = ['cxcalc', '-S', '--ignore-error', 'majormicrospecies', '-H', '7.4', '-K',
               f'{"-M" if tautomerize else ""}', input_fname]
    with open(output_fname, 'w') as file:
        subprocess.run(cmd_run, stdout=file, text=True)


def read_protonate_chemaxon(fname):
    for mol in Chem.SDMolSupplier(fname, sanitize=False):
        if mol:
            mol_name = mol.GetProp('_Name')
            smi = mol.GetPropsAsDict().get('MAJORMS', None)
            if smi is not None:
                yield smi, mol_name



def protonate_dimorphite(input_fname, output_fname):
    executable = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dimorphite_dl', 'dimorphite_dl.py')
    dimorphite_run(smiles_file=input_fname, output_file=output_fname, max_variants=1, silent=True, min_ph=7.4, max_ph=7.4)


def read_smiles(fname):
    with open(fname) as f:
        for line in f:
            yield tuple(line.strip().split()[:2])