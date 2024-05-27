import os
import re
import sys
import traceback
from multiprocessing import cpu_count

from meeko import (MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy,
                   RDKitMolCreate)
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from rdkit.Chem.rdMolAlign import AlignMolConformers
from typing import Optional, Iterator

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
    pdbqt_string_list = []
    preparator = MoleculePreparation(hydrate=False, flexible_amides=False, rigid_macrocycles=True, min_ring_size=7, max_ring_size=33)
    try:
        cids = [x.GetId() for x in mol.GetConformers()]
        for cid in cids:
            mol_setups = preparator.prepare(mol, conformer_id=cid)

            for setup in mol_setups:
                pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
                if not is_ok:
                    print(f"{mol.GetProp('_Name')} has error in converting to pdbqt: {error_msg}")
                    return None
                
                pdbqt_string_list.append(pdbqt_string)
                if verbose:
                    print(f"{setup}")
    except Exception:
        sys.stderr.write(
            "Warning. Incorrect mol object to convert to pdbqt. Continue. \n"
        )
        traceback.print_exc()
        pdbqt_string_list = None

    return pdbqt_string_list


def GetConformerRMSFromAtomIds(mol: Chem.Mol, confId1: int, confId2: int, atomIds: Optional[list[int]]=None, prealigned: bool=False) -> float:
    """ Returns the RMS between two conformations based on the atomIds passed as the input.
        By default, the conformers will be aligned to the first conformer
        before the RMS calculation and, as a side-effect, the second will be left
        in the aligned state.

    Arguments:
        - mol:        the molecule
        - confId1:    the id of the first conformer
        - confId2:    the id of the second conformer
        - atomIds:    (optional) list of atom ids to use a points for
                        alingment **AND RMSD** calculation- defaults to all atoms
        - prealigned: (optional) by default the conformers are assumed
                        be unaligned and the second conformer be aligned
                        to the first

    """
  # align the conformers if necessary
  # Note: the reference conformer is always the first one
    if not prealigned:
        if atomIds:
            AlignMolConformers(mol, confIds=[confId1, confId2], atomIds=atomIds)
        else:
            AlignMolConformers(mol, confIds=[confId1, confId2])

  # calculate the RMS between the two conformations
    conf1 = mol.GetConformer(id=confId1)
    conf2 = mol.GetConformer(id=confId2)
    ssr = 0
    if atomIds:
        for i in atomIds:
            d = conf1.GetAtomPosition(i).Distance(conf2.GetAtomPosition(i))
            ssr += d * d
        ssr /= len(atomIds)
    else:
        for i in range(mol.GetNumAtoms()):
            d = conf1.GetAtomPosition(i).Distance(conf2.GetAtomPosition(i))
            ssr += d * d
        ssr /= mol.GetNumAtoms()

    return np.sqrt(ssr)

def GetConformerRMSMatrixForSaturatedRingMolecule(mol: Chem.Mol, atomIds:list[list[int]]=None, prealigned: bool=False) -> list[float]:
    """ Returns the RMS matrix of the conformers of a molecule based on the alignment and rmsd of saturated ring.
        The function calculates the mean of the RMSD of each saturated ring with the GetConformerRMSFromAtomIds.
        The alignment is done per ring (for example, three alignments are done for three saturated ring) 
        
        By default, the conformers will be aligned to the first conformer
        before the RMS calculation and, as a side-effect, the second will be left
        in the aligned state.

    As a side-effect, the conformers will be aligned to the first
    conformer (i.e. the reference) and will left in the aligned state.

    Arguments:
        - mol:     the molecule
        - atomIds: (optional) list of atom ids to use a points for
                    alingment - defaults to all atoms
        - prealigned: (optional) by default the conformers are assumed
                        be unaligned and will therefore be aligned to the
                        first conformer

    Note that the returned RMS matrix is symmetrical, i.e. it is the
    lower half of the matrix, e.g. for 5 conformers::

        rmsmatrix = [ a,
                      b, c,
                      d, e, f,
                      g, h, i, j]

    where a is the RMS between conformers 0 and 1, b is the RMS between
    conformers 0 and 2, etc.
    This way it can be directly used as distance matrix in e.g. Butina
    clustering.

    """

    cmat_list = []
    total_ring_atoms = sum([len(atom_in_ring) for atom_in_ring in atomIds])
    for atom_id_in_ring in atomIds:
        # if necessary, align the conformers
        # Note: the reference conformer is always the first one
        rmsvals = []
        confIds = [conf.GetId() for conf in mol.GetConformers()]
        if not prealigned:
            if atom_id_in_ring:
                AlignMolConformers(mol, atomIds=atom_id_in_ring, RMSlist=rmsvals)
            else:
                AlignMolConformers(mol, RMSlist=rmsvals)
        else:  # already prealigned
            for i in range(1, len(confIds)):
                rmsvals.append(GetConformerRMSFromAtomIds(mol, confIds[0], confIds[i], atomIds=atom_id_in_ring, prealigned=prealigned))
        # loop over the conformations (except the reference one)
        cmat_per_ring = []
        for i in range(1, len(confIds)):
            cmat_per_ring.append(rmsvals[i - 1])
            for j in range(1, i):
                cmat_per_ring.append(GetConformerRMSFromAtomIds(mol, confIds[i], confIds[j], atomIds=atom_id_in_ring, prealigned=True))

        cmat_list.append(np.square(np.array(cmat_per_ring))*len(atom_id_in_ring))

    cmat_list_array = np.array(cmat_list)
    return list(np.sqrt(np.sum(cmat_list_array, axis=0)/total_ring_atoms))
  

def mol_embedding_3d(mol: Chem.Mol, seed: int=43) -> Chem.Mol:

    def find_saturated_ring_with_substituent(mol: Chem.Mol) -> list[list[int]]:
        ssr = Chem.GetSymmSSSR(mol)
        saturated_ring_list = []
        for ring_idx in ssr:
            is_atom_saturated_array = [mol.GetAtomWithIdx(atom_id).GetHybridization() == Chem.HybridizationType.SP3 for atom_id in ring_idx]
            if any(is_atom_saturated_array):
                ring_and_substituent_idx = []
                for ring_atom_idx in ring_idx:
                    ring_and_substituent_idx += [x.GetIdx() for x in mol.GetAtomWithIdx(ring_atom_idx).GetNeighbors()]
                ring_and_substituent_idx = list(set(ring_and_substituent_idx))
                saturated_ring_list.append(ring_and_substituent_idx)

        #to show that the function pick the correct neighbors and inspect visually.
        #to be removed before merging the PR
        img = Draw.MolsToGridImage([mol],highlightAtomLists=[list(set([ x for ring in saturated_ring_list for x in ring]))])
        img.save(f'images/{mol}_{mol.GetProp("_Name")}.png')
        return saturated_ring_list

    def gen_conf(mole: Chem.Mol, useRandomCoords: bool, randomSeed: int, has_saturated_ring: bool) -> tuple[Chem.Mol, float]:
        params = AllChem.ETKDGv3()
        params.useRandomCoords = useRandomCoords
        params.randomSeed = randomSeed
        if has_saturated_ring:
            #10 is used as default according to the C language documentation iirc, but I have to specify the numbers.
            conf_stat = AllChem.EmbedMultipleConfs(mole, 100, params)
        else:
            conf_stat = AllChem.EmbedMolecule(mole, params)
        return mole, conf_stat
    
    def remove_confs_rms(mol: Chem.Mol, saturated_ring_list: list[list[int]], rms: float=0.25, keep_nconf: Optional[int]=None) -> Chem.Mol:
        """
        The function uses AgglomerativeClustering to select conformers.

        :param mol: input molecule with multiple conformers
        :param rms: discard conformers which are closer than given value to a kept conformer
        :param keep_nconf: keep at most the given number of conformers. This parameter has precedence over rms
        :return:
        """

        def gen_ids(ids: int) -> Iterator[int]:
            for i in range(1, len(ids)):
                for j in range(0, i):
                    yield j, i

        if keep_nconf and mol.GetNumConformers() <= keep_nconf:
            return mol

        if mol.GetNumConformers() <= 1:
            return mol

        rms_ = GetConformerRMSMatrixForSaturatedRingMolecule(Chem.RemoveHs(mol), atomIds=saturated_ring_list, prealigned=False)

        cids = [c.GetId() for c in mol.GetConformers()]
        arr = np.zeros((len(cids), len(cids)))
        for (i, j), v in zip(gen_ids(cids), rms_):
            arr[i, j] = v
            arr[j, i] = v

        cl = AgglomerativeClustering(n_clusters=None, linkage='complete', metric='precomputed', distance_threshold=rms).fit(arr)
        keep_ids = []
        for i in set(cl.labels_):
            ids = np.where(cl.labels_ == i)[0]
            j = arr[np.ix_(ids, ids)].mean(axis=0).argmin()
            keep_ids.append(cids[ids[j]])
        remove_ids = set(cids) - set(keep_ids)

        for cid in sorted(remove_ids, reverse=True):
            mol.RemoveConformer(cid)

        while True:
            ids = np.in1d(cids, keep_ids)
            arr = arr[np.ix_(ids, ids)]

            #sometimes clustering result gives matrix < rms when rms is high enough
            if all(arr[arr != 0] < rms) or not any(arr[arr != 0] < rms):
                break

            cids = [c.GetId() for c in mol.GetConformers()]

            #keep ids with maximum rms
            max_rms_ids = np.unravel_index(np.argmax(arr, axis=None), arr.shape)
            count_below_rms = np.count_nonzero((arr < rms) & (arr !=0), axis=1)

            max_count_below_rms_ids = np.argwhere(count_below_rms == np.amax(count_below_rms))
            remove_ids = np.setdiff1d(max_count_below_rms_ids, np.array(max_rms_ids))

            # remove max rms if remove id is empty
            if len(remove_ids) ==0:
                remove_ids = np.setdiff1d(max_count_below_rms_ids, np.array([]))
            remove_cids = [cids[id] for id in remove_ids]

            keep_ids = np.array(list(set(cids) - set(list(remove_cids))))

            for cid in sorted(remove_cids, reverse=True):
                mol.RemoveConformer(cid)

        if keep_nconf and mol.GetNumConformers() > keep_nconf and mol.GetNumConformers() > 1:
            ids = np.in1d(cids, keep_ids)
            arr = arr[np.ix_(ids, ids)]   # here other indexing operation should be used, because ids is a boolean array
            cl = AgglomerativeClustering(n_clusters=keep_nconf, linkage='complete', metric='precomputed').fit(arr)

            keep_ids = []
            cids = [c.GetId() for c in mol.GetConformers()]
            for i in set(cl.labels_):
                ids = np.where(cl.labels_ == i)[0]
                j = arr[np.ix_(ids, ids)].mean(axis=0).argmin()
                keep_ids.append(cids[ids[j]])
            remove_ids = set(cids) - set(keep_ids)

            for cid in sorted(remove_ids, reverse=True):
                mol.RemoveConformer(cid)

        return mol
        
    if not isinstance(mol, Chem.Mol):
        return None
    
    saturated_rings_with_substituents = find_saturated_ring_with_substituent(mol)
    has_saturated_ring = (len(saturated_rings_with_substituents)>0)

    mol = Chem.AddHs(mol, addCoords=True)
    if not mol_is_3d(mol):  # only for non 3D input structures
        mol, conf_stat = gen_conf(mol, useRandomCoords=False, randomSeed=seed, has_saturated_ring=has_saturated_ring)
        if conf_stat == -1:
            # if molecule is big enough and rdkit cannot generate a conformation - use params.useRandomCoords = True
            mol, conf_stat = gen_conf(mol, useRandomCoords=True, randomSeed=seed, has_saturated_ring=has_saturated_ring)
            if conf_stat == -1:
                return None
        AllChem.UFFOptimizeMolecule(mol, maxIters=100)

        writer = Chem.SDWriter(f'conformer/{mol.GetProp("_Name")}_before_remove.sdf')
        for cid in [c.GetId() for c in mol.GetConformers()]:
            writer.write(mol, confId=cid)
        print(f"[For Testing Only] {mol.GetProp('_Name')} has {len(saturated_rings_with_substituents)} saturated ring")
        print(f"[For Testing Only] Before removing conformation: {mol.GetProp('_Name')} has {mol.GetNumConformers()} conf")
        mol = remove_confs_rms(mol, saturated_rings_with_substituents, rms=0.5)
        print(f"[For Testing Only] After removing conformation: {mol.GetProp('_Name')} has {mol.GetNumConformers()} conf")
        
        AlignMolConformers(mol)

        writer = Chem.SDWriter(f'conformer3/{mol.GetProp("_Name")}_after_remove_050.sdf')
        for cid in [c.GetId() for c in mol.GetConformers()]:
            writer.write(mol, confId=cid)
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
        sys.stderr.write(f"Parsing PDB was failed: {mol_id}\n")
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
