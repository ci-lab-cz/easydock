import subprocess
import sys

from functools import partial
from multiprocessing import Pool
from rdkit import Chem
from dimorphite_dl.dimorphite_dl import run as dimorphite_run
from math import ceil
import tempfile
import os

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

def chunk_into_n(smi_l: list[str], n: int):
    smi_size = ceil(len(smi_l) / n)
    return list(map(lambda x: smi_l[x * smi_size:x * smi_size + smi_size], list(range(n))))

def dummy_protonate_dimorphite(input_output_fname: tuple[str, str]):
    input_fname, output_fname = input_output_fname
    dimorphite_run(smiles_file=input_fname, output_file=output_fname, max_variants=1, silent=True, min_ph=7.4, max_ph=7.4)

def protonate_dimorphite(input_fname: str, output_fname: str, ncpu: int = 1):

    with open(input_fname,'r') as input_file:
        smi_l = input_file.readlines()

    chunk_smi_l = chunk_into_n(smi_l, n=ncpu)

    temp_fname_list = []
    for chunk in chunk_smi_l:
        with tempfile.NamedTemporaryFile(suffix='.smi',mode='w',encoding='utf-8', delete=False) as input_tmp:

            input_tmp.write(''.join(chunk))
            output_tmp = tempfile.NamedTemporaryFile(suffix='.smi',mode='w',encoding='utf-8', delete=False)
            temp_fname_list.append((input_tmp.name, output_tmp.name))
    
    pool = Pool(ncpu)
    pool.map(dummy_protonate_dimorphite, temp_fname_list)
    pool.close()
    pool.join()

    with open(output_fname, 'wt') as output_file:
        for temp_fname in temp_fname_list:
            input_temp_fname, output_temp_fname = temp_fname

            with open(output_temp_fname,'r') as output_smi:
                output_file.write(''.join(output_smi))
            
            os.remove(input_temp_fname)
            os.remove(output_temp_fname)

def read_smiles(fname):
    with open(fname) as f:
        for line in f:
            yield tuple(line.strip().split()[:2])
