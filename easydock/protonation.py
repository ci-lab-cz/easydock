import contextlib
import os
import subprocess
import tempfile

from functools import partial
from math import ceil
from multiprocessing import Pool
from typing import Iterator

from rdkit import Chem
from easydock.auxiliary import chunk_into_n
from easydock.read_input import read_input

"""
There two types of protonation programs:

1. File-based approaches. They require an input file and return an output file. The integration should include two 
implemented functions: 
1) protonate_xxx with two necessary arguments input_fname and output_fname (input and output filenames). 
Input file will be SMILES with names separated by tab. Output file format may be any (output is a temporary file which 
does not have an extension).
2) read_protonate_xxx takes a file name with protonated molecules in the format corresponding to 
the protonate_xxx function and returns a generator of tuples (SMILES, mol_name).

2. Python-based approaches. They utilize pure Python workflow and should use multiprocessing.pool to make enumeration 
efficient. The implementation consists of a single function:
1) protonate_xxx which is a generator. It should take items argument which is a generator over tuples of (smi, mol_name) 
and yield a tuple of (SMILES, mol_name).

These functions should be intergated in database.add_protonation function, there is a special section of initialization 
of protonation functions. All functions may take additional arguments, which should be passed with 
partial(protonate_xxx, arg1=value1, ...) at the intialization step.   
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


def __protonate_dimorphite_mp(input_output_fname: tuple[str, str]):
    from dimorphite_dl.dimorphite_dl import run as dimorphite_run
    input_fname, output_fname = input_output_fname
    dimorphite_run(smiles_file=input_fname, output_file=output_fname, max_variants=1, silent=True, min_ph=7.4, max_ph=7.4)


def protonate_dimorphite(input_fname: Iterator[tuple[str, str]], output_fname: str, ncpu: int = 1):

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
    pool.map(__protonate_dimorphite_mp, temp_fname_list)
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


def protonate_pkasolver(items: str, ncpu: int = 1, smi_size=1):
    import torch
    from pkasolver.query import QueryModel

    chunksize = min(max(1, smi_size // ncpu), 500)
    model = QueryModel()
    with contextlib.redirect_stdout(None):
        if torch.cuda.is_available() or ncpu == 1:
            for smi, mol_name in items:
                yield __protonate_pkasolver((smi, mol_name), model=model)
        else:
            pool = Pool(ncpu)
            for smi, mol_name in pool.imap_unordered(partial(__protonate_pkasolver, model=model), items, chunksize=chunksize):
                yield smi, mol_name


def __protonate_pkasolver(args, model):
    from pkasolver.query import calculate_microstate_pka_values
    smi, mol_name = args
    mol = Chem.MolFromSmiles(smi, sanitize=True)
    ph = 7.4
    states = calculate_microstate_pka_values(mol, only_dimorphite=False, query_model=model)
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
        ids = output_mol.GetSubstructMatches(Chem.MolFromSmarts('[$([NH+,NH2+,NH3+]-[*]=O)]'))
        if ids:
            for i in ids:
                output_mol.GetAtomWithIdx(i[0]).SetFormalCharge(0)

    return Chem.MolToSmiles(output_mol), mol_name


def protonate_molgpka(items: str, ncpu: int = 1, smi_size=1):
    from molgpka.predict_pka_mp import load_state_dicts

    models = load_state_dicts()
    chunksize = min(max(1, smi_size // ncpu), 500)

    pool = Pool(ncpu)
    for smi, mol_name in pool.imap_unordered(partial(__protonate_molgpka, models=models), items, chunksize=chunksize):
        yield smi, mol_name


def add_hydrogen_to_atom(editable_mol, atom_idx):
    hydrogen = Chem.Atom(1)
    hydrogen_idx = editable_mol.AddAtom(hydrogen)
    editable_mol.AddBond(atom_idx, hydrogen_idx, Chem.BondType.SINGLE)


def __protonate_molgpka(args, models):
    from molgpka.predict_pka_mp import predict
    from rdkit.Chem import AllChem
    from copy import deepcopy

    smi, mol_name = args
    #print(smi, mol_name)
    mol = Chem.MolFromSmiles(smi, sanitize=True)
    ph = 7.4
    try:
        #editable_mol = Chem.RWMol(mol)
        #print(Chem.MolToSmiles(editable_mol), mol_name)
        base_dict, acid_dict, mol_ = predict(mol, models, uncharged=True)
        editable_mol = Chem.RWMol(mol_)
        for idx, pKa in sorted(acid_dict.items(), reverse=True):
            atom = editable_mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 1:
                atom_idx = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()][0]
                atom_not_H = editable_mol.GetAtomWithIdx(atom_idx)
                charge = atom_not_H.GetFormalCharge()
                Hs = atom_not_H.GetTotalNumHs(includeNeighbors=True)
                if pKa < ph:
                    if Hs > 0 and charge >= 0:
                        editable_mol.RemoveAtom(idx)
                        atom_not_H.SetFormalCharge(charge - 1)
                        atom_not_H.UpdatePropertyCache()
                    else:
                        print('mol with problem to deprotonate', mol_name)

        for idx, pKa in base_dict.items():
            atom = editable_mol.GetAtomWithIdx(idx)
            charge = atom.GetFormalCharge()
            if pKa > ph:
                if charge <= 0:
                    add_hydrogen_to_atom(editable_mol, idx)
                    atom.SetFormalCharge(charge + 1)
                else:
                    print('mol with problem to protonate', mol_name)
                atom.UpdatePropertyCache()

        changed_smi = Chem.MolToSmiles(Chem.RemoveHs(editable_mol))
        # print(f'{datetime.datetime.now()}; PID: {os.getpid()}; {item[1]} end')
        # sys.stdout.flush()
    except Exception as e:
        print(e)
        print(mol_name)
        return None, mol_name
    return changed_smi, mol_name
