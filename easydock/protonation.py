import contextlib
import logging
import os
import subprocess
import tempfile
import traceback
import warnings

from functools import partial
from multiprocessing import Pool
from typing import Iterator, Tuple

from rdkit import Chem
from easydock.containers import apptainer_exec
from easydock.auxiliary import chunk_into_n, expand_path

"""
There two types of protonation programs:

1. File-based approaches. They require an input file and return an output file. The integration should include two 
implemented functions: 
1) protonate_xxx with two necessary arguments input_fname and output_fname (input and output filenames). 
Input file will be SMILES with names separated by tab. Output file format may be any (output is a temporary file which 
does not have an extension).
2) read_protonate_xxx takes a file name with protonated molecules in the format corresponding to 
the protonate_xxx function and returns a generator of tuples (SMILES, mol_name).

There is a special case of file-based approaches if they use a specific python environment. In this case they should be 
implemented inside an apptainer container (.sif file), where all necessary dependencies will be installed. 
The container should have a command "protonate" which takes two necessary arguments -i/--input and -o/--output, 
where input and output files can be passed. Thus, the protonation can be invoked by a command:
apptainer run container.sif protonate -i input.smi -o output.smi
Integration of these approaches into easydock is similar as described above. The difference is that protonate_xxx 
function should take an additional argument - path to apptainer container.

2. Python-based approaches. They utilize pure Python workflow and should use multiprocessing.pool to make enumeration 
efficient. The implementation consists of a single function:
1) protonate_xxx which is a generator. It should take items argument which is a generator over tuples of (smi, mol_name) 
and yield a tuple of (SMILES, mol_name).

These functions should be integrated in database.add_protonation function, there is a special section of initialization 
of protonation functions. All functions may take additional arguments, which should be passed with 
partial(protonate_xxx, arg1=value1, ...) at the intialization step.   
"""

# MolGpKa fix patterns
# patterns + ids of atoms in a pattern to be fixed + type of a pattern
# this group of patterns has two centers and only one with more favorable pKa/pKb value will be remained after fix
molgpka_patterns1 = [('[#7+;!$([#7+][O-]);!$([#7+H0v4])]~*~[#7+;!$([#7+][O-]);!$([#7+H0v4])]', (0, 2), 'pos'),    # 21
                     ('[#7+;!$([#7+][O-]);!$([#7+H0v4])]**[#7+;!$([#7+][O-]);!$([#7+H0v4])]', (0, 3), 'pos'),   # 3
                     ('[#7+;!$([#7+][O-]);!$([#7+H0v4])]***[#7+;!$([#7+][O-]);!$([#7+H0v4])]', (0, 4), 'pos'),  # 4
                     ('c[n-][n-]c', (1, 2), 'neg'),     # 9
                     ('[O-]cc[O-]', (0, 3), 'neg'),     # 13
                     ('[O-]cac[O-]', (0, 4), 'neg'),    # 14
                     ('[O-]caac[O-]', (0, 5), 'neg')]   # 15
molgpka_patterns1 = [(Chem.MolFromSmarts(p), ids, p_type) for p, ids, p_type in molgpka_patterns1]

# this group of patterns has one center which will be reverted to a neutral state
molgpka_patterns2 = [('[$([NH+,NH2+,NH3+]-[*]=[O,S])]', 'pos'), # 1
                     ('[$([N-](C(=O))[N-]S(=O)=O)]', 'neg'),  # 2
                     ('[$([N-]1C(=O)N(C)C(=O)C1)]', 'neg'),  # 5
                     ('[$([NH-](C(=O))C1=NC=CS1)]', 'neg'),  # 6
                     ('[$([NH2+]C(=[NH+])[NH2])]', 'pos'),  # 7
                     ('[$([N-]C(=[NH+])[NH2])]', 'neg'),  # 8
                     ('[$([nH+]1aaaa1)]', 'pos'),  # 10
                     ('[$([n-]1c(=O)[n-]c(=O)cn1)]', 'neg'),  # 11
                     ('[$([n-]1ncnc1[N-]S(=O)=O)]', 'neg'),  # 12
                     ('[$([O-]aaC(=O)[O-])]', 'neg'),  # 16
                     ('[$([#8-][#6][#7-])]', 'neg'),  # 17
                     ('[$([O-][CX4][CX3][O-])]', 'neg'),  # 18
                     ('[$([N+]=C[N-]S(=O)=O)]', 'pos'),  # 19
                     ('[$([n-]1cnc(=O)c2cccn12)]', 'neg'),  # 20
                     ('[$([N-](a)C([NH2])=[NH2+])]', 'neg'),  #22
                     ('[$([NX4+;H1,H2,H3]c)]', 'pos'), #23
                     ('[$([O-][CX4])]', 'neg')]  #24
molgpka_patterns2 = [(Chem.MolFromSmarts(p), p_type) for p, p_type in molgpka_patterns2]


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


def __protonate_dimorphite_mp(input_output_fname: Tuple[str, str]):
    from dimorphite_dl.dimorphite_dl import run as dimorphite_run
    input_fname, output_fname = input_output_fname
    dimorphite_run(smiles_file=input_fname, output_file=output_fname, max_variants=1, silent=True, min_ph=7.4, max_ph=7.4)


def protonate_dimorphite(input_fname: Iterator[Tuple[str, str]], output_fname: str, ncpu: int = 1):

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


def protonate_pkasolver(items: Iterator[Tuple[str, str]], ncpu: int = 1, mol_count=1):
    import torch
    from pkasolver.query import QueryModel

    chunksize = min(max(1, mol_count // ncpu), 500)
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


def protonate_molgpka(items: Iterator[Tuple[str, str]], ncpu: int = 1):
    # parallel execution of protonation was disabled because runs much slower than a single process protonation
    warnings.filterwarnings('ignore', category=UserWarning)
    from molgpka.predict_pka_mp import load_state_dicts, load_models
    warnings.filterwarnings("default", category=UserWarning)

    models = load_state_dicts()
    models = load_models(models)
    for q in items:
        yield __protonate_molgpka(q, models)


def __add_hydrogen_to_atom(editable_mol, atom_idx):
    hydrogen = Chem.Atom(1)
    hydrogen_idx = editable_mol.AddAtom(hydrogen)
    editable_mol.AddBond(atom_idx, hydrogen_idx, Chem.BondType.SINGLE)


def __assign_pka_pkb_to_heavy_atoms(mol, acid_dict, base_dict):
    for idx, value in acid_dict.items():
        atom_idx = [neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors()][0]
        mol.GetAtomWithIdx(atom_idx).SetDoubleProp('pka', float(value))
    for idx, value in base_dict.items():
        mol.GetAtomWithIdx(idx).SetDoubleProp('pkb', float(value))
    return mol


def __protonate_molgpka(args, models):
    from molgpka.predict_pka_mp import predict2

    smi, mol_name = args
    mol = Chem.MolFromSmiles(smi, sanitize=True)
    ph = 7.4

    warnings.filterwarnings('ignore', category=UserWarning)

    try:

        # mol_ will have al H explicit (acidic pKa is assigned to H)
        base_dict, acid_dict, mol_ = predict2(mol, models, uncharged=True)
        mol_ = __assign_pka_pkb_to_heavy_atoms(mol_, acid_dict, base_dict)
        editable_mol = Chem.RWMol(Chem.RemoveHs(mol_))

        # assign protonation states
        for atom in editable_mol.GetAtoms():
            if atom.HasProp('pka') and atom.GetDoubleProp('pka') < ph:
                Hs = atom.GetTotalNumHs()
                charge = atom.GetFormalCharge()
                if Hs > 0 and charge >= 0:
                    atom.SetFormalCharge(charge - 1)
                    # after RemoveHs such hydrogens remain explicit, thus we have to reduce their number explicitly
                    # however, if there is no explicit Hs, implicit Hs should be recalculated with UpdatePropertyCache
                    if atom.GetNumExplicitHs() > 0:
                        atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
                    atom.UpdatePropertyCache()
                else:
                    logging.warning(f'(molgpka) Molecule {mol_name} has issues with assignment of protonation states (deprotonation)')
            elif atom.HasProp('pkb') and atom.GetDoubleProp('pkb') > ph:
                charge = atom.GetFormalCharge()
                if charge <= 0:
                    atom.SetFormalCharge(charge + 1)
                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                    atom.UpdatePropertyCache()
                else:
                    logging.warning(f'(molgpka) Molecule {mol_name} has issues with assignment of protonation states (protonation)')
        editable_mol.UpdatePropertyCache()

        # fix and revert some protonation states
        for pattern, ids, pattern_type in molgpka_patterns1:
            matches = editable_mol.GetSubstructMatches(pattern, uniquify=True)
            if pattern_type == 'pos':
                prop_name = 'pkb'
            else:
                prop_name = 'pka'
            while matches:
                worst_value = float('inf') if pattern_type == 'pos' else float('-inf')
                worst_atom_idx = None
                for match in matches:
                    v1 = editable_mol.GetAtomWithIdx(match[ids[0]]).GetDoubleProp(prop_name)
                    v2 = editable_mol.GetAtomWithIdx(match[ids[1]]).GetDoubleProp(prop_name)
                    if pattern_type == 'neg' and max(v1, v2) > worst_value:
                        worst_value = max(v1, v2)
                        if v1 > v2:
                            worst_atom_idx = match[ids[0]]
                        else:
                            worst_atom_idx = match[ids[1]]
                    elif pattern_type == 'pos' and min(v1, v2) < worst_value:
                        worst_value = min(v1, v2)
                        if v1 < v2:
                            worst_atom_idx = match[ids[0]]
                        else:
                            worst_atom_idx = match[ids[1]]

                charge = editable_mol.GetAtomWithIdx(worst_atom_idx).GetFormalCharge()
                if pattern_type == 'pos':
                    atom = editable_mol.GetAtomWithIdx(worst_atom_idx)
                    atom.SetFormalCharge(charge - 1)
                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
                elif pattern_type == 'neg':
                    atom = editable_mol.GetAtomWithIdx(worst_atom_idx)
                    atom.SetFormalCharge(charge + 1)
                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)

                matches = editable_mol.GetSubstructMatches(pattern, uniquify=True)

        editable_mol.UpdatePropertyCache()

        # single atom fixes
        for pattern, pattern_type in molgpka_patterns2:
            matches = editable_mol.GetSubstructMatches(pattern, uniquify=True)
            for match in matches:
                atom = editable_mol.GetAtomWithIdx(match[0])
                if pattern_type == 'pos':
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
                    atom.UpdatePropertyCache()
                if pattern_type == 'neg':
                    atom.SetFormalCharge(atom.GetFormalCharge() + 1)
                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                    atom.UpdatePropertyCache()
        editable_mol.UpdatePropertyCache()

        changed_smi = Chem.MolToSmiles(Chem.RemoveHs(editable_mol))

    except Exception as e:
        logging.warning(f'{mol_name} caused an error during protonation, it will be skipped\n'
                        f'{traceback.format_exc()}')
        return None, mol_name

    finally:
        warnings.filterwarnings("default", category=UserWarning)

    return changed_smi, mol_name


def protonate_apptainer(input_fname: str, output_fname: str, container_fname: str) -> None:

    bind_path = set()
    bind_path.add(os.path.dirname(expand_path(input_fname)))
    bind_path.add(os.path.dirname(expand_path(output_fname)))

    apptainer_exec(container_fname,
                   ['protonate', '-i', input_fname, '-o', output_fname],
                   list(bind_path))
