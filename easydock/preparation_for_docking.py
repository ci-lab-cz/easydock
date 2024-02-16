import os
import re
import sys
import traceback
from multiprocessing import cpu_count

from meeko import (MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy,
                   RDKitMolCreate)
from rdkit import Chem
from rdkit.Chem import AllChem


def cpu_type(x):
    return max(1, min(int(x), cpu_count()))


def filepath_type(x):
    if x:
        return os.path.abspath(x)
    else:
        return x


def mol_is_3d(mol):
    if mol.GetConformers() and list(mol.GetConformers())[0].Is3D():
        return True
    return False


def mol_from_smi_or_molblock(ligand_string):
    mol = Chem.MolFromMolBlock(ligand_string)
    if mol is None:
        mol = Chem.MolFromSmiles(ligand_string)
    return mol


def mk_prepare_ligand(mol, verbose=False):
    preparator = MoleculePreparation(hydrate=False, flexible_amides=False, rigid_macrocycles=True, min_ring_size=7, max_ring_size=33)
    try:
        mol_setups = preparator.prepare(mol)
        for setup in mol_setups:
            pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
            if verbose:
                print(f"{setup}")
    except Exception:
        sys.stderr.write(
            "Warning. Incorrect mol object to convert to pdbqt. Continue. \n"
        )
        traceback.print_exc()
        pdbqt_string = None

    return pdbqt_string


def mol_embedding_3d(mol, seed=43):

    def gen_conf(mole, useRandomCoords, randomSeed):
        params = AllChem.ETKDGv3()
        params.useRandomCoords = useRandomCoords
        params.randomSeed = randomSeed
        conf_stat = AllChem.EmbedMolecule(mole, params)
        return mole, conf_stat

    if not isinstance(mol, Chem.Mol):
        return None
    mol = Chem.AddHs(mol, addCoords=True)
    if not mol_is_3d(mol):  # only for non 3D input structures
        mol, conf_stat = gen_conf(mol, useRandomCoords=False, randomSeed=seed)
        if conf_stat == -1:
            # if molecule is big enough and rdkit cannot generate a conformation - use params.useRandomCoords = True
            mol, conf_stat = gen_conf(mol, useRandomCoords=True, randomSeed=seed)
            if conf_stat == -1:
                return None
        AllChem.UFFOptimizeMolecule(mol, maxIters=100)
    return mol


def ligand_preparation(mol, boron_replacement=False, seed=43):
    """
    If input ligand is not a 3D structure a conformer will be generated by RDKit, otherwise the provided 3D structure
    will be used. Boron atoms are replaced with carbon to enable docking using Vina and gnina
    :param mol:
    :param boron_replacement: indicate to whether replace boron with carbon atoms or not
    :param seed: fixed to 43 to generate consistent random stereoisomers for compounds with undefined stereocenters
    :return: PDBQT block
    """

    try:
        mol = mol_embedding_3d(mol, seed=seed)
        if mol:
            if boron_replacement:
                idx_boron = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 5]
                for id_ in idx_boron:
                    if mol.GetAtomWithIdx(id_).GetFormalCharge() < 0:
                        mol.GetAtomWithIdx(id_).SetFormalCharge(0)
                    mol.GetAtomWithIdx(id_).SetAtomicNum(6)
                    # mol.UpdatePropertyCache() # uncomment if necessary
            mol_conf_pdbqt = mk_prepare_ligand(mol, verbose=False)
            return mol_conf_pdbqt
    except Exception:
        traceback.print_exc()
        return None


def assign_bonds_from_template(template_mol, mol):
    # explicit hydrogends are removed from carbon atoms (chiral hydrogens) to match pdbqt mol,
    # e.g. [NH3+][C@H](C)C(=O)[O-]
    template_mol_ = Chem.Mol(template_mol)
    template_mol_ = Chem.AddHs(template_mol_, explicitOnly=True,
                               onlyOnAtoms=[a.GetIdx() for a in template_mol_.GetAtoms() if
                                            a.GetAtomicNum() != 6])
    mol = AllChem.AssignBondOrdersFromTemplate(template_mol_, mol)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    return mol


def boron_reduction(mol_B, mol):
    mol_B_ = Chem.Mol(mol_B)
    mol_ = Chem.Mol(mol)

    idx_boron = {a.GetIdx(): a.GetFormalCharge() for a in mol_B_.GetAtoms() if a.GetAtomicNum() == 5}
    if idx_boron:

        for id_, charge in idx_boron.items():
            if charge < 0:
                mol_B_.GetAtomWithIdx(id_).SetFormalCharge(0)
            mol_B_.GetAtomWithIdx(id_).SetAtomicNum(6)

        mol_ = assign_bonds_from_template(mol_B_, mol_)
        idx = mol_.GetSubstructMatches(mol_B_)
        mol_idx_boron = [tuple(sorted((ids[i], j) for i, j in idx_boron.items())) for ids in idx]
        mol_idx_boron = list(set(mol_idx_boron))  # retrieve all ids matched possible boron atom positions
        if len(mol_idx_boron) == 1:  # check whether this set of ids is unique
            for id_, charge in mol_idx_boron[0]:
                mol_.GetAtomWithIdx(id_).SetAtomicNum(5)
                mol_.GetAtomWithIdx(id_).SetFormalCharge(charge)
        else:  # if not - several equivalent mappings exist
            sys.stderr.write('different mappings was detected. The structure cannot be recostructed automatically.')
            return None

    return mol_


def pdbqt2molblock(pdbqt_block, template_mol, mol_id):
    """
    The function takes PDBQT block with one or more poses and converts top pose to MDL MOL format. The function tries
    to return back boron atoms
    :param pdbqt_block: a single string with a single PDBQT block (a single pose)
    :param template_mol: Mol of a reference structure to assign bond orders
    :param mol_id: name of a molecule which will be added as a title in the output MOL block
    :param boron_replacement: indicate whether to try to return boron atoms instead af carbon ones
    :return: a single string with a MOL block, if conversion failed returns None
    """
    mol_block = None
    try:
        pdbqt_mol = PDBQTMolecule(pdbqt_block, is_dlg=False, skip_typing=True, poses_to_read=1)
        rdkitmol_list = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        rdkit_mol = rdkitmol_list[0]
    except Exception:
        traceback.print_exc()
        sys.stderr.write(f"Parsing PDB was failed (fixing did not help): {mol_id}\n")
        return None

    try:
        if 5 in [atom.GetAtomicNum() for atom in template_mol.GetAtoms()]:
            mol = boron_reduction(template_mol, rdkit_mol)
        else:
            mol = assign_bonds_from_template(template_mol, rdkit_mol)

        mol.SetProp("_Name", mol_id)
        mol_block = Chem.MolToMolBlock(mol)
    except Exception:
        traceback.print_exc()
        sys.stderr.write(f"Could not assign bond orders while parsing PDB: {mol_id}. Trying to fix.\n")

    return mol_block
