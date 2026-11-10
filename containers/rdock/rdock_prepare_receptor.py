#!/usr/bin/env python3
"""
Receptor preparation for rDock: a protonated protein PDB -> MOL2 receptor, system definition (.prm) and cavity (.as).

Used by rdock_server.py (automatic preparation at `init`) and as a command line tool
(`apptainer run rdock.sif prepare ...` or `python3 rdock_prepare_receptor.py ...`).

Steps:
  1. ATOM/HETATM records of the first model are kept (waters dropped unless --keep-waters) and written to a clean PDB.
  2. Open Babel converts the PDB to MOL2 (bond perception, Tripos atom types).
  3. The MOL2 is corrected, because rDock derives its receptor model from the file and Open Babel gets three things
     wrong on real proteins:
       - residue identity: subst_id becomes the residue position in the file and a SUBSTRUCTURE section carries the
         chain as the rDock segment (Open Babel writes the residue number, which merges residues of different chains
         and of insertion codes into one group, and rDock then assigns wrong group charges);
       - aromatic side-chain rings of HIS, TRP, PHE and TYR are typed C.ar/N.ar by name (Open Babel perceives
         aromaticity from geometry and misses a part of them, which gives rDock phantom implicit hydrogens and wrong
         acceptors);
       - every BOND record is written with the lower atom id first (rDock's flexible-receptor code picks the side
         of an OH/NH3+ bond to rotate from the atom order, and on a receptor made of several fragments the wrong
         order rotates a whole fragment).
     Coordinates, atom names and partial charges are not touched. Atom count and order must match the PDB.
  4. A .prm with the receptor and the cavity mapper (reference ligand or a sphere) is written, rbcavity computes the
     cavity. rbcavity exits 0 on some fatal errors, so its log decides.
  5. Optional pharmacophore restraints (SECTION PHARMA, a penalty on poses that miss the given points): the points
     are rewritten beside the .prm as pharma_mandatory.const / pharma_optional.const, so the receptor folder does not
     depend on the original restraint files.

The protein must already carry hydrogens with standard PDB atom names (PDBFixer, reduce, Open Babel is not enough):
rDock takes polar hydrogens as donors and assigns group charges by residue and atom names (LYS HZ1-3, ...).
"""
import argparse
import logging
import os
import re
import shutil
import subprocess
import sys

RING_TYPES = {
    'HIS': {'CG': 'C.ar', 'CD2': 'C.ar', 'CE1': 'C.ar', 'ND1': 'N.ar', 'NE2': 'N.ar'},
    'TRP': {**dict.fromkeys(('CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'), 'C.ar'), 'NE1': 'N.ar'},
    'PHE': dict.fromkeys(('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'), 'C.ar'),
}
RING_TYPES['TYR'] = RING_TYPES['PHE']

WATERS = {'HOH', 'WAT', 'DOD', 'H2O', 'TIP', 'TIP3', 'SOL'}

# rbcavity/rbdock print a fatal error as "NAME at <file>, line N" at the start of a line and may still exit 0.
# A harmless warning starts with "**WARNING**", so it does not match "^\w+ at".
BAD_LOG = re.compile(r'Mismatch between SUBSTRUCTURE|Inconsistent subunit names|Removing orphan hydrogen|'
                     r'^\w+ at \S+, line \d+(?:\n[^\n]*)?', re.M)

MOL2_NAME = 'protein.mol2'
REFLIG_NAME = 'ref_ligand.sdf'

# pharmacophore restraint types of rDock (RbtConstraint.cxx): Don matches the donor hydrogen, Aro a ring centre
PHARMA_TYPES = ('Any', 'Don', 'Acc', 'Aro', 'Hyd', 'Hal', 'Har', 'Ani', 'Cat')
PHARMA_FILES = {'constraints_file': 'pharma_mandatory.const', 'optional_file': 'pharma_optional.const'}
PHARMA_KEYS = ('constraints_file', 'optional_file', 'nopt', 'weight', 'write_errors')

logger = logging.getLogger(__name__)


def read_pdb_atoms(path, keep_waters=False):
    """ATOM/HETATM records of the first model in file order:
    list of (line, atom_name, res_name, chain, res_seq, ins_code)."""
    atoms = []
    with open(path) as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                res_name = line[17:20].strip()
                if not keep_waters and res_name in WATERS:
                    continue
                atoms.append((line.rstrip('\n'), line[12:16].strip(), res_name,
                              line[21].strip(), line[22:26].strip(), line[26].strip()))
            elif line.startswith('ENDMDL'):
                break
    if not atoms:
        raise ValueError(f'no ATOM/HETATM records in {path}')
    return atoms


def has_hydrogens(atoms):
    for line, name, *_ in atoms:
        element = line[76:78].strip() if len(line) >= 78 else ''
        if element == 'H' or (not element and name[:1] == 'H'):
            return True
    return False


def write_pdb(atoms, path):
    with open(path, 'w') as f:
        for atom in atoms:
            f.write(atom[0] + '\n')
        f.write('END\n')


def mol2_fixup(atoms, raw_mol2, out_mol2):
    """Open Babel MOL2 of the PDB atoms -> rDock receptor MOL2 (see the module docstring). Returns the residue count."""
    positions, previous, n = [], None, 0
    for atom in atoms:
        key = atom[2:6]  # res_name, chain, res_seq, ins_code
        if key != previous:
            n += 1
            previous = key
        positions.append(n)

    # sections by header; only MOLECULE, ATOM and BOND are kept (Open Babel may add UNITY_ATOM_ATTR with formal
    # charges from the PDB charge column; rDock assigns charges from its own residue table and ignores them)
    sections, current = {}, None
    with open(raw_mol2) as f:
        for ln in f.read().split('\n'):
            if ln.startswith('@<TRIPOS>'):
                current = ln.strip()
                sections[current] = []
            elif current is not None:
                sections[current].append(ln)
    try:
        mol_lines = sections['@<TRIPOS>MOLECULE']
        atom_lines = [ln for ln in sections['@<TRIPOS>ATOM'] if ln.strip()]
        bond_lines = [ln for ln in sections['@<TRIPOS>BOND'] if ln.strip()]
    except KeyError as e:
        raise ValueError(f'{raw_mol2} has no {e} section')
    if len(atom_lines) != len(atoms):
        raise ValueError(f'Open Babel wrote {len(atom_lines)} atoms for {len(atoms)} PDB atoms')

    out_atoms, residues = [], []
    for k, ln in enumerate(atom_lines):
        t = ln.split()
        _, pdb_name, res_name, chain, res_seq, ins_code = atoms[k]
        if t[1] != pdb_name:
            raise ValueError(f'atom {k + 1}: MOL2 name {t[1]} differs from PDB name {pdb_name} (atom order changed)')
        sub_name = f'{res_name}{res_seq}{ins_code}'
        if positions[k] > len(residues):
            residues.append((positions[k], sub_name, int(t[0]), chain or '_', res_name))
        atom_type = RING_TYPES.get(res_name, {}).get(pdb_name, t[5])
        out_atoms.append(f'{int(t[0]):7d} {t[1]:<4s} {t[2]:>12s}{t[3]:>10s}{t[4]:>10s} '
                         f'{atom_type:<5s}{positions[k]:5d}  {sub_name:<8s}{t[8]:>10s}')

    bonds = []
    for ln in bond_lines:
        bond_id, p, q, order = ln.split()[:4]
        p, q = int(p), int(q)
        bonds.append(f'{int(bond_id):6d}{min(p, q):6d}{max(p, q):6d} {order}')

    counts = mol_lines[1].split()  # title, counts, mol type, charge type, ...
    mol_lines[1] = f' {counts[0]} {counts[1]} {len(residues)} 0 0'
    substructures = [f'{i:6d} {name:<8s}{root:7d} RESIDUE           1 {chain} {res_name}     0'
                     for i, name, root, chain, res_name in residues]
    with open(out_mol2, 'w') as f:
        f.write('\n'.join(['@<TRIPOS>MOLECULE'] + mol_lines + ['@<TRIPOS>ATOM'] + out_atoms + ['@<TRIPOS>BOND']
                          + bonds + ['@<TRIPOS>SUBSTRUCTURE'] + substructures) + '\n')
    return len(residues)


def mapper_section(reflig_name=None, center=None, radius=None, small_sphere=None, large_sphere=4.0,
                   min_volume=100, max_cavities=1, gridstep=0.5):
    """Lines of the SECTION MAPPER block: the reference ligand method (reflig_name) or the two-sphere method (center).
    Defaults follow rDock's examples: radius 6.0 and small sphere 1.0 around a ligand, radius 10.0 and small sphere 1.5
    around a point."""
    if reflig_name:  # the ligand mapper has no LARGE_SPHERE parameter (rbcavity: RBT_BAD_ARGUMENT)
        lines = ['    SITE_MAPPER RbtLigandSiteMapper',
                 f'    REF_MOL {reflig_name}',
                 f'    RADIUS {6.0 if radius is None else float(radius)}',
                 f'    SMALL_SPHERE {1.0 if small_sphere is None else float(small_sphere)}']
    elif center is not None:
        x, y, z = (float(v) for v in center)
        lines = ['    SITE_MAPPER RbtSphereSiteMapper',
                 f'    CENTER ({x:.3f},{y:.3f},{z:.3f})',
                 f'    RADIUS {10.0 if radius is None else float(radius)}',
                 f'    SMALL_SPHERE {1.5 if small_sphere is None else float(small_sphere)}',
                 f'    LARGE_SPHERE {float(large_sphere)}']
    else:
        raise ValueError('either a reference ligand or a center is required to define the cavity')
    lines += [f'    MIN_VOLUME {int(min_volume)}',
              f'    MAX_CAVITIES {int(max_cavities)}',
              f'    GRIDSTEP {float(gridstep)}']
    return lines


def write_prm(path, mapper_lines, receptor_flex=None, title='easydock', pharma_lines=None):
    """System definition file. File names inside are relative to its folder: rDock splits values on whitespace, so an
    absolute path with a space would break, and rbcavity/rbdock are run with that folder as the working directory."""
    lines = ['RBT_PARAMETER_FILE_V1.00', f'TITLE {title}', f'RECEPTOR_FILE {MOL2_NAME}']
    if receptor_flex:
        lines.append(f'RECEPTOR_FLEX {float(receptor_flex)}')
    lines += ['SECTION MAPPER'] + mapper_lines + ['END_SECTION',
                                                   'SECTION CAVITY',
                                                   '    SCORING_FUNCTION RbtCavityGridSF',
                                                   '    RMAX 0.1',
                                                   'END_SECTION']
    if pharma_lines:
        lines += pharma_lines
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def parse_restraints(value):
    """Pharmacophore restraint points as a list of (x, y, z, tolerance, type) from a restraints file path, from the
    text of such a file (several lines) or from a list of [x, y, z, tolerance, type] items ('x y z tolerance type'
    strings are accepted in the list too). As in rDock, fields are separated by whitespace or commas; blank lines and
    '#' comments are skipped (a MOE .ph4 file is therefore not accepted)."""
    if isinstance(value, str):
        path = os.path.expanduser(value)
        if os.path.isfile(path):
            with open(path) as f:
                value = f.read()
        elif '\n' not in value.strip() and len(value.split()) != 5:  # one line of five fields is a restraint
            raise ValueError(f'restraints file not found: {value}')
        items = value.splitlines()
    elif isinstance(value, (list, tuple)):
        items = list(value)
    else:
        raise ValueError(f'restraints must be a file path or a list of [x, y, z, tolerance, type], '
                         f'got {type(value).__name__}')
    points = []
    for item in items:
        if isinstance(item, str):
            item = item.split('#', 1)[0].replace(',', ' ').split()
            if not item:
                continue
        if not isinstance(item, (list, tuple)) or len(item) != 5:
            raise ValueError(f'restraint {item!r}: expected x y z tolerance type')
        try:
            x, y, z, tol = (float(v) for v in item[:4])
        except (TypeError, ValueError):
            raise ValueError(f'restraint {item!r}: x, y, z and the tolerance must be numbers')
        kind = str(item[4])
        if kind not in PHARMA_TYPES:
            raise ValueError(f'restraint {item!r}: the type must be one of {", ".join(PHARMA_TYPES)} (case matters)')
        if tol <= 0:
            raise ValueError(f'restraint {item!r}: the tolerance must be positive')
        points.append((x, y, z, tol, kind))
    if not points:
        raise ValueError('no restraint points found')
    return points


def write_restraints(points, path):
    with open(path, 'w') as f:
        f.write(''.join(f'{x:.3f} {y:.3f} {z:.3f} {tol:.2f} {kind}\n' for x, y, z, tol, kind in points))


def pharma_section(constraints_name, optional_name=None, nopt=0, weight=1.0, write_errors=False):
    """Lines of the SECTION PHARMA block (rDock reference guide, 'Pharmacophore restraints')."""
    lines = ['SECTION PHARMA',
             '    SCORING_FUNCTION RbtPharmaSF',
             f'    WEIGHT {float(weight)}',
             f'    CONSTRAINTS_FILE {constraints_name}']
    if optional_name:
        lines += [f'    OPTIONAL_FILE {optional_name}', f'    NOPT {int(nopt)}']
    if write_errors:
        lines.append('    WRITE_ERRORS TRUE')
    return lines + ['END_SECTION']


def prepare_pharma(pharma, out_dir):
    """Writes the restraint files of the `pharma` mapping (keys constraints_file, optional_file, nopt, weight,
    write_errors; the two files in any form parse_restraints accepts) into out_dir under fixed names and returns the
    SECTION PHARMA lines referring to them."""
    if not isinstance(pharma, dict):
        raise ValueError("'pharma' must be a mapping with constraints_file (required), optional_file, nopt, weight, "
                         "write_errors")
    unknown = sorted(set(pharma) - set(PHARMA_KEYS))
    if unknown:
        raise ValueError(f'unknown pharma keys {unknown}; allowed: {", ".join(PHARMA_KEYS)}')
    if pharma.get('constraints_file') is None:
        raise ValueError('pharma.constraints_file (the mandatory restraints) is required')
    names, counts = {}, {}
    for key, fname in PHARMA_FILES.items():
        if pharma.get(key) is not None:
            points = parse_restraints(pharma[key])
            write_restraints(points, os.path.join(out_dir, fname))
            names[key], counts[key] = fname, len(points)
    nopt = int(pharma.get('nopt') or 0)
    n_optional = counts.get('optional_file', 0)
    if not 0 <= nopt <= n_optional:
        raise ValueError(f'pharma.nopt must be between 0 and the number of optional restraints ({n_optional}), '
                         f'got {nopt}')
    if n_optional and nopt == 0:
        logger.warning('pharma.nopt is 0: the %d optional restraints have no effect', n_optional)
    weight = float(pharma.get('weight', 1.0))
    if weight < 0:
        raise ValueError('pharma.weight must be >= 0')
    return pharma_section(names['constraints_file'], names.get('optional_file'), nopt=nopt, weight=weight,
                          write_errors=bool(pharma.get('write_errors', False)))


def write_reference_ligand(src, dst):
    """The first record of the reference ligand file with Kekulé bond orders: rDock reads bond orders 1-3 only, and an
    aromatic bond type (4, as Marvin writes it) gives 'makes a bond with an unknown bond order' errors. Without RDKit,
    or when RDKit cannot read the file, it is copied as it is."""
    try:
        from rdkit import Chem
    except ImportError:
        shutil.copyfile(src, dst)
        return
    mol = next((m for m in Chem.SDMolSupplier(src, removeHs=False) if m is not None), None)
    if mol is None:
        shutil.copyfile(src, dst)
        return
    Chem.MolToMolFile(mol, dst)


def run_logged(cmd, log_path, cwd, timeout):
    with open(log_path, 'w') as log:
        proc = subprocess.run(cmd, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, text=True, timeout=timeout)
    with open(log_path, errors='replace') as f:
        return proc.returncode, f.read()


def prepare(protein, out_dir, reflig=None, center=None, name='receptor', radius=None, small_sphere=None,
            large_sphere=4.0, min_volume=100, max_cavities=1, gridstep=0.5, receptor_flex=None, keep_waters=False,
            allow_no_hydrogens=False, pharma=None, timeout=1800):
    """Writes out_dir/{protein.mol2, <name>.prm, <name>.as, ref_ligand.sdf} and the logs, plus the restraint files
    when `pharma` (a mapping, see prepare_pharma) is given. Returns a dict with the prm path, the cavity volume (A^3,
    or None if rbcavity did not print it) and the residue count.
    Raises ValueError/RuntimeError with a message naming the failed step."""
    for tool in ('obabel', 'rbcavity'):
        if shutil.which(tool) is None:
            raise RuntimeError(f'{tool} not found on PATH')
    if '.' in name:
        raise ValueError(f'name {name!r} must not contain a dot: rDock derives the cavity file name from the .prm '
                         f'name up to its first dot')
    protein = os.path.abspath(os.path.expanduser(protein))
    if not os.path.isfile(protein):
        raise ValueError(f'protein PDB not found: {protein}')
    out_dir = os.path.abspath(os.path.expanduser(out_dir))
    os.makedirs(out_dir, exist_ok=True)

    atoms = read_pdb_atoms(protein, keep_waters=keep_waters)
    if not allow_no_hydrogens and not has_hydrogens(atoms):
        raise ValueError(f'{protein} has no hydrogens; rDock needs a protonated receptor with standard atom names '
                         f'(PDBFixer, reduce). Use allow_no_hydrogens to override.')
    clean_pdb = os.path.join(out_dir, f'{name}_input.pdb')
    write_pdb(atoms, clean_pdb)

    raw_mol2 = os.path.join(out_dir, f'{name}_obabel.mol2')
    rc, log = run_logged(['obabel', '-ipdb', clean_pdb, '-omol2', '-O', raw_mol2],
                         os.path.join(out_dir, f'{name}_obabel.log'), out_dir, timeout)
    if rc or not os.path.isfile(raw_mol2) or not os.path.getsize(raw_mol2):
        raise RuntimeError(f'obabel failed (exit {rc}), see {name}_obabel.log in {out_dir}')
    n_residues = mol2_fixup(atoms, raw_mol2, os.path.join(out_dir, MOL2_NAME))
    os.unlink(raw_mol2)

    reflig_name = None
    if reflig:
        reflig = os.path.abspath(os.path.expanduser(reflig))
        if not os.path.isfile(reflig):
            raise ValueError(f'reference ligand not found: {reflig}')
        reflig_name = REFLIG_NAME
        write_reference_ligand(reflig, os.path.join(out_dir, REFLIG_NAME))
    mapper = mapper_section(reflig_name=reflig_name, center=center, radius=radius, small_sphere=small_sphere,
                            large_sphere=large_sphere, min_volume=min_volume, max_cavities=max_cavities,
                            gridstep=gridstep)
    pharma_lines = prepare_pharma(pharma, out_dir) if pharma else None
    prm = os.path.join(out_dir, f'{name}.prm')
    write_prm(prm, mapper, receptor_flex=receptor_flex, pharma_lines=pharma_lines)

    cavity = os.path.join(out_dir, f'{name}.as')
    if os.path.exists(cavity):
        os.unlink(cavity)  # a stale cavity must not pass for this run's
    rc, log = run_logged(['rbcavity', '-r', f'{name}.prm', '-W'], os.path.join(out_dir, f'{name}_rbcavity.log'),
                         out_dir, timeout)
    bad = BAD_LOG.search(log)
    if rc or bad or not os.path.isfile(cavity) or not os.path.getsize(cavity):
        reason = f'log has {bad.group(0).replace(chr(10), ": ")!r}' if bad else f'no {name}.as written'
        raise RuntimeError(f'rbcavity failed (exit {rc}): {reason}, see {name}_rbcavity.log in {out_dir}')
    volume = re.search(r'Total volume\s+([\d.]+)', log)
    return {'prm': prm,
            'cavity_volume': float(volume.group(1)) if volume else None,
            'n_residues': n_residues}


def main(argv=None):
    parser = argparse.ArgumentParser(description='Prepare an rDock receptor (MOL2, .prm, cavity) from a protonated '
                                                 'protein PDB. The cavity is defined by a reference ligand or a '
                                                 'sphere around a point.')
    parser.add_argument('-p', '--protein', required=True, help='protein PDB with hydrogens and standard atom names')
    site = parser.add_mutually_exclusive_group(required=True)
    site.add_argument('-l', '--reflig', help='reference ligand SDF (3D) defining the cavity')
    site.add_argument('--center', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                      help='center of a sphere defining the cavity')
    parser.add_argument('-o', '--out_dir', required=True, help='output folder')
    parser.add_argument('--name', default='receptor', help='stem of the .prm and .as files (default: receptor)')
    parser.add_argument('--radius', type=float,
                        help='mapping radius around the reference ligand (default 6.0) or the center (default 10.0)')
    parser.add_argument('--small_sphere', type=float, help='small probe radius (default 1.0 ligand, 1.5 sphere)')
    parser.add_argument('--large_sphere', type=float, default=4.0,
                        help='large probe radius of the sphere mapper (default 4.0; not used with --reflig)')
    parser.add_argument('--min_volume', type=int, default=100, help='minimum cavity volume in A^3 (default 100)')
    parser.add_argument('--max_cavities', type=int, default=1, help='number of cavities to keep (default 1)')
    parser.add_argument('--gridstep', type=float, default=0.5, help='mapping grid step (default 0.5)')
    parser.add_argument('--receptor_flex', type=float,
                        help='make receptor OH/NH3+ groups within this distance of the cavity flexible (e.g. 3.0); '
                             'default: rigid receptor')
    parser.add_argument('--keep_waters', action='store_true', help='keep water molecules in the receptor')
    parser.add_argument('--allow_no_hydrogens', action='store_true',
                        help='do not stop when the protein has no hydrogens (not recommended)')
    parser.add_argument('--pharma_constraints', metavar='FILE',
                        help='rDock pharmacophore restraints file with the mandatory points, one "x y z tolerance '
                             f'type" per line, types {", ".join(PHARMA_TYPES)}; adds SECTION PHARMA to the .prm')
    parser.add_argument('--pharma_optional', metavar='FILE', help='restraints file with the optional points')
    parser.add_argument('--pharma_nopt', type=int, default=0,
                        help='number of optional points a pose has to meet (default 0)')
    parser.add_argument('--pharma_weight', type=float, default=1.0,
                        help='weight of the restraint penalty in the total score (default 1.0)')
    parser.add_argument('--pharma_write_errors', action='store_true',
                        help='rbdock writes ligands lacking the mandatory features to <output>_errors.sd')
    args = parser.parse_args(argv)
    pharma = None
    if args.pharma_constraints:
        pharma = {'constraints_file': args.pharma_constraints, 'optional_file': args.pharma_optional,
                  'nopt': args.pharma_nopt, 'weight': args.pharma_weight, 'write_errors': args.pharma_write_errors}
    elif args.pharma_optional:
        parser.error('--pharma_optional requires --pharma_constraints')
    try:
        res = prepare(args.protein, args.out_dir, reflig=args.reflig, center=args.center, name=args.name,
                      radius=args.radius, small_sphere=args.small_sphere, large_sphere=args.large_sphere,
                      min_volume=args.min_volume, max_cavities=args.max_cavities, gridstep=args.gridstep,
                      receptor_flex=args.receptor_flex, keep_waters=args.keep_waters,
                      allow_no_hydrogens=args.allow_no_hydrogens, pharma=pharma)
    except (ValueError, RuntimeError, subprocess.TimeoutExpired) as e:
        sys.exit(f'rdock receptor preparation failed: {e}')
    volume = f'{res["cavity_volume"]} A^3' if res['cavity_volume'] is not None else 'unknown'
    print(f'receptor: {res["prm"]} ({res["n_residues"]} residues, cavity volume {volume})')


if __name__ == '__main__':
    main()
