import contextlib
import os
import subprocess
import sys
import tempfile

from math import ceil
from multiprocessing import Pool

from rdkit import Chem
from read_input import read_input

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


def protonate_dimorphite_mp(input_output_fname: tuple[str, str]):
    from dimorphite_dl.dimorphite_dl import run as dimorphite_run
    input_fname, output_fname = input_output_fname
    dimorphite_run(smiles_file=input_fname, output_file=output_fname, max_variants=1, silent=True, min_ph=7.4, max_ph=7.4)


def protonate_dimorphite(input_fname: str, output_fname: str, ncpu: int = 1):

    with open(input_fname,'r') as input_file:
        smi_l = input_file.readlines()

    chunk_smi_l = chunk_into_n(smi_l, n=ncpu)

    temp_fname_list = []
    for chunk in chunk_smi_l:
        with tempfile.NamedTemporaryFile(suffix='.smi', mode='w', encoding='utf-8', delete=False) as input_tmp:
            input_tmp.write(''.join(chunk))
            output_tmp = tempfile.NamedTemporaryFile(suffix='.smi', mode='w', encoding='utf-8', delete=False)
            temp_fname_list.append((input_tmp.name, output_tmp.name))
    
    pool = Pool(ncpu)
    pool.map(protonate_dimorphite_mp, temp_fname_list)
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


class DummyFile(object):
    def write(self, x): pass


@contextlib.contextmanager
def nostd():
    save_stdout = sys.stdout
    save_stderr = sys.stderr
    sys.stdout = DummyFile()
    sys.stderr = DummyFile()
    yield
    sys.stdout = save_stdout
    sys.stderr = save_stderr


def protonate_pkasolver(input_fname: str, output_fname: str, ncpu: int = 1):
    from pkasolver.query import QueryModel
    model = QueryModel()
    pool = Pool(ncpu)
    with open(output_fname, 'wt') as f:
        for smi, name in pool.imap_unordered(__protonate_pkasolver,
                                             ((mol, mol_name, model) for (mol, mol_name) in read_input(input_fname))):
            f.write(f'{smi}\t{name}\n')


def __protonate_pkasolver(args):
    from pkasolver.query import calculate_microstate_pka_values
    mol, mol_name, model = args
    ph = 7.4
    with nostd():
        states = calculate_microstate_pka_values(mol, only_dimorphite=True, query_model=model)
    if not states:
        output_mol = mol
    else:
        # select protonated form of the first state with pKa > pH or deprotonated form of the state with the highest pKa
        output_mol = None
        for state in states:
            if state.pka > ph:
                output_mol = state.protonated_mol
                break
        if not output_mol:
            output_mol = states[-1].deprotonated_mol
        # fix protonated amides and analogs (about 2% such structures in 23K molecules)
        output_mol = Chem.RemoveHs(output_mol)
        ids = output_mol.GetSubstructMatches(Chem.MolFromSmarts('[$([NH+]-[*]=O)]'))
        if ids:
            for i in ids:
                output_mol.GetAtomWithIdx(i[0]).SetFormalCharge(0)
    return Chem.MolToSmiles(output_mol), mol_name
