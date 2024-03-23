import os.path
import subprocess
import sys

from functools import partial
from multiprocessing import Pool
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

def chunk_into_n_file(fname, number_of_files):
    with open(fname) as infp:
        fname_list = [fname[:-4]+'_input'+str(i)+'.smi' for i in range(number_of_files)]
        files = [open(fname[:-4]+'_input'+str(i)+'.smi','w') for i in range(number_of_files)]
        for i, line in enumerate(infp):
            files[i % number_of_files].write(line)
        for f in files:
            f.close()
        return fname_list

def dummy_protonate_dimorphite(input_fname):
    dimorphite_run(smiles_file=input_fname, output_file=input_fname[:-4]+'_output.smi', max_variants=1, silent=True, min_ph=7.4, max_ph=7.4)

def protonate_dimorphite(input_fname, output_fname, ncpu):
    
    pool = Pool(ncpu)
    fname_list = chunk_into_n_file(input_fname, ncpu)
    
    with open(output_fname, 'wt') as output_file:
        pool.map(dummy_protonate_dimorphite, fname_list)

        for fname in fname_list:
            output_file.writelines(open(fname[:-4]+'_output.smi','r').readlines())

def read_smiles(fname):
    with open(fname) as f:
        for line in f:
            yield tuple(line.strip().split()[:2])
