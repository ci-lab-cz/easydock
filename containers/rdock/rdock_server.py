#!/usr/bin/env python3
"""
rDock docking server for EasyDock (`--program server`), JSON Lines over stdin/stdout (docs/server_protocol.md).

  info  formats: ligands arrive as 3D mol blocks (`mol3d`: EasyDock embeds them with hydrogens), poses leave as SDF,
        lower score is better, one molecule per request by default
  init  the receptor, either prepared here from a protonated PDB and a reference ligand or a center (see
        rdock_prepare_receptor.py) into a temporary folder of this worker, optionally with pharmacophore restraints
        (`pharma`), or a ready system definition (.prm) whose cavity (.as) lies beside it; plus the docking protocol
        and settings
  dock  one SD file and one rbdock run per request; per ligand the poses are sorted by `sort_field` (SCORE) and the
        first pose's `score_field` (SCORE.INTER) is the docking score, as rDock's authors recommend: the total score
        picks the pose, the intermolecular score ranks ligands

The formal charges of a ligand are copied from the `M  CHG` lines into the atom block, the only place rDock reads them.
rbdock exits 0 on fatal errors, so the request fails on an error header in its output without any pose written.
A ligand lacking the features of the mandatory pharmacophore restraints is rejected by rbdock (RBT_LIGAND_ERROR):
it is returned as failed and the reason is logged, the request does not fail.
rDock docks only the largest fragment of a ligand and drops the hydrogens bonded to carbon; the returned pose gets
those hydrogens back with RDKit when it is available.
"""
import atexit
import json
import logging
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile

try:
    from rdkit import Chem
except ImportError:  # the pose is then returned as rbdock wrote it
    Chem = None

import rdock_prepare_receptor as receptor_prep

logger = logging.getLogger(__name__)

INFO = {
    # one molecule per rbdock run: rbdock seeds its random generator once per run, so a larger batch makes a
    # ligand's poses depend on the batch it was in. Raise batch_size via info_server for throughput.
    'batch_size': 1,
    'ligand_in_format': 'mol3d',
    'ligand_out_format': 'sdf',
    'score_mode': 'min',
}

# MDL atom block charge codes
CHG = {3: 1, 2: 2, 1: 3, 0: 0, -1: 5, -2: 6, -3: 7}
# a fatal rbdock error: "NAME at <file>, line N" at the start of a line, usually followed by a message line
RBDOCK_ERROR = re.compile(r'^(\w+) at \S+, line \d+(?:\n([^\n]*))?', re.M)


def _expand(path):
    if path is None:
        return None
    return os.path.abspath(os.path.expanduser(str(path)))


def ligand_record(mol_block, title):
    """A mol block from EasyDock -> one SD record for rbdock: the title set, the coordinates centred (rDock keeps the
    input position and splits the atom line on whitespace, so a crystal frame far from the origin can run two columns
    together) and the formal charges from `M  CHG` written into the atom block."""
    lines = mol_block.rstrip('\n').split('\n')
    if len(lines) < 4:
        raise ValueError('not a mol block')
    counts = lines[3]
    if 'V3000' in counts:
        raise ValueError('V3000 mol block: rDock reads V2000 only')
    n_atoms, n_bonds = int(counts[0:3]), int(counts[3:6])
    atom_lines = [ln.ljust(48) for ln in lines[4:4 + n_atoms]]
    if len(atom_lines) != n_atoms:
        raise ValueError('truncated mol block')
    coords = [(float(ln[0:10]), float(ln[10:20]), float(ln[20:30])) for ln in atom_lines]
    if lines[1][20:22] != '3D' and all(abs(z) < 1e-6 for _, _, z in coords):
        raise ValueError('2D ligand: rDock needs 3D coordinates (the server expects ligand_in_format mol3d)')
    charges = {}
    for ln in lines[4 + n_atoms + n_bonds:]:
        if ln.startswith('M  CHG'):
            fields = ln.split()
            for i in range(int(fields[2])):
                charges[int(fields[3 + 2 * i])] = int(fields[4 + 2 * i])
    centre = [sum(c[i] for c in coords) / n_atoms for i in range(3)]
    out = []
    for i, ln in enumerate(atom_lines):
        charge = charges.get(i + 1, 0)
        if charge not in CHG:
            raise ValueError(f'formal charge {charge} cannot be written in the MDL atom block')
        x, y, z = (coords[i][k] - centre[k] for k in range(3))
        out.append(f'{x:10.4f}{y:10.4f}{z:10.4f}{ln[30:36]}{CHG[charge]:3d}{ln[39:]}'.rstrip())
    lines[0] = title
    return '\n'.join(lines[:4] + out + lines[4 + n_atoms:]) + '\n$$$$\n'


def split_sd(text):
    """SD records (each ending with the $$$$ line) in file order."""
    records = []
    for chunk in text.split('$$$$\n'):
        if chunk.strip():
            records.append(chunk + '$$$$\n')
    return records


def sd_title(record):
    return record.split('\n', 1)[0].strip()


def sd_field(record, name):
    m = re.search(r'^>\s*<' + re.escape(name) + r'>[^\n]*\n([^\n]*)', record, re.M)
    return m.group(1).strip() if m else None


def pose_mol_block(record, title):
    """The pose as a mol block titled `title`: the record up to M  END, with the hydrogens rDock dropped added back
    when RDKit is available."""
    block = record.split('M  END', 1)[0] + 'M  END\n'
    block = title + block[block.find('\n'):]
    if Chem is None:
        return block
    try:
        mol = Chem.MolFromMolBlock(block, removeHs=False)
        if mol is None:
            return block
        mol = Chem.AddHs(mol, addCoords=True)
        Chem.AssignStereochemistryFrom3D(mol)
        mol.SetProp('_Name', title)
        return Chem.MolToMolBlock(mol)
    except Exception:
        logger.warning('%s: could not add hydrogens to the pose, returning it as rbdock wrote it', title)
        return block


class _RDockServer:
    def __init__(self):
        self.initialized = False
        self.prm_dir = None
        self.prm_name = None
        self.receptor_tmp = None
        self.protocol = 'dock.prm'
        self.n_runs = 50
        self.n_poses = 5
        self.seed = 0
        self.sort_field = 'SCORE'
        self.score_field = 'SCORE.INTER'
        self.extra_args = []
        self.tmpdir = None
        self.keep_tmpdir = False
        atexit.register(self.cleanup)

    def cleanup(self):
        if self.receptor_tmp and not self.keep_tmpdir:
            shutil.rmtree(self.receptor_tmp, ignore_errors=True)
        self.receptor_tmp = None

    def _resolve_protocol(self, protocol):
        path = _expand(protocol)
        if os.path.isfile(path):
            return path
        rbt_root = os.environ.get('RBT_ROOT')
        if rbt_root and not os.path.isfile(os.path.join(rbt_root, 'data', 'scripts', protocol)):
            raise ValueError(f'protocol {protocol!r} is neither a file nor in {rbt_root}/data/scripts')
        return protocol  # a name rbdock looks up in $RBT_ROOT/data/scripts (dock.prm, dock_solv.prm, ...)

    def init(self, payload):
        try:
            return self._init(payload)
        except Exception as e:
            logger.exception('init failed')
            return {'status': 'error', 'error': f'{type(e).__name__}: {e}'}

    def _init(self, payload):
        self.initialized = False
        self.cleanup()
        if shutil.which('rbdock') is None:
            return {'status': 'error', 'error': 'rbdock not found on PATH'}
        p = dict(payload or {})
        self.tmpdir = _expand(p.get('tmpdir'))
        self.keep_tmpdir = bool(p.get('keep_tmpdir', False))
        self.protocol = self._resolve_protocol(str(p.get('protocol', 'dock.prm')))
        self.n_runs = int(p.get('n_runs', 50))
        self.n_poses = int(p.get('n_poses', 5))
        seed = p.get('seed', 0)
        self.seed = None if seed is None else int(seed)  # seed: null lets rbdock pick a random seed
        self.sort_field = str(p.get('sort_field', 'SCORE'))
        self.score_field = str(p.get('score_field', 'SCORE.INTER'))
        extra = p.get('extra_args', [])
        self.extra_args = shlex.split(extra) if isinstance(extra, str) else [str(a) for a in extra]
        if self.n_runs < 1 or self.n_poses < 1:
            return {'status': 'error', 'error': 'n_runs and n_poses must be positive'}

        receptor = p.get('receptor')
        protein = p.get('protein')
        pharma = p.get('pharma')
        if receptor and protein:
            return {'status': 'error', 'error': "give either 'receptor' (a ready .prm) or 'protein' (a PDB to prepare)"}
        if pharma is not None and not isinstance(pharma, dict):
            return {'status': 'error', 'error': "'pharma' must be a mapping with constraints_file (required), "
                                                "optional_file, nopt, weight, write_errors"}
        if receptor:
            prm = _expand(receptor)
            if not os.path.isfile(prm):
                return {'status': 'error', 'error': f'receptor .prm not found: {prm}'}
            self.prm_dir, self.prm_name = os.path.split(prm)
            cavity = os.path.join(self.prm_dir, self.prm_name.split('.')[0] + '.as')
            if not os.path.isfile(cavity) or not os.path.getsize(cavity):
                return {'status': 'error',
                        'error': f'cavity file {cavity} not found: run `rbcavity -r {self.prm_name} -W` in '
                                 f'{self.prm_dir} or prepare the receptor with rdock_prepare_receptor.py'}
            if p.get('receptor_flex') is not None:
                logger.warning("'receptor_flex' is ignored with a ready .prm: set RECEPTOR_FLEX in the file")
            if pharma:
                return {'status': 'error', 'error': f"'pharma' cannot be combined with a ready .prm: add SECTION PHARMA "
                                                    f"to {prm} (rDock reference guide, 'Pharmacophore restraints') or "
                                                    f"give 'protein' to prepare the receptor here"}
        elif protein:
            self.receptor_tmp = tempfile.mkdtemp(prefix='rdock_receptor_', dir=self.tmpdir)
            res = receptor_prep.prepare(
                protein=_expand(protein), out_dir=self.receptor_tmp, reflig=_expand(p.get('reflig')),
                center=p.get('center'), name='receptor', radius=p.get('radius'), small_sphere=p.get('small_sphere'),
                large_sphere=p.get('large_sphere', 4.0), min_volume=p.get('min_volume', 100),
                max_cavities=p.get('max_cavities', 1), gridstep=p.get('gridstep', 0.5),
                receptor_flex=p.get('receptor_flex'), keep_waters=bool(p.get('keep_waters', False)),
                allow_no_hydrogens=bool(p.get('allow_no_hydrogens', False)), pharma=pharma or None)
            self.prm_dir, self.prm_name = self.receptor_tmp, 'receptor.prm'
            logger.info('receptor prepared in %s: %d residues, cavity volume %s A^3',
                        self.receptor_tmp, res['n_residues'], res['cavity_volume'])
        else:
            return {'status': 'error', 'error': "'receptor' (a ready .prm) or 'protein' (a PDB) is required"}

        self.initialized = True
        return {'status': 'ok'}

    def dock(self, payload):
        if not self.initialized:
            return {'status': 'error', 'error': 'Server not initialized'}
        if not isinstance(payload, dict):
            return {'status': 'error', 'error': f'payload must be a dict {{name: mol_block}}, got {type(payload)}'}
        results = {name: None for name in payload}
        workdir = tempfile.mkdtemp(prefix='rdock_dock_', dir=self.tmpdir)
        try:
            records, sent = [], []
            for name, block in payload.items():
                try:
                    if not isinstance(block, str) or not block.strip():
                        raise ValueError('empty ligand')
                    records.append(ligand_record(block, name))
                    sent.append(name)
                except Exception as e:
                    logger.warning('%s skipped: %s', name, e)
            if not records:
                return {'status': 'ok', 'results': results}

            sd_in = os.path.join(workdir, 'ligands.sd')
            root = os.path.join(workdir, 'docked')
            with open(sd_in, 'w') as f:
                f.write(''.join(records))
            cmd = ['rbdock', '-i', sd_in, '-o', root, '-r', self.prm_name, '-p', self.protocol, '-n', str(self.n_runs)]
            if self.seed is not None:
                cmd += ['-s', str(self.seed)]
            cmd += self.extra_args
            logger.debug('Running in %s: %s', self.prm_dir, ' '.join(cmd))
            proc = subprocess.run(cmd, cwd=self.prm_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            with open(os.path.join(workdir, 'rbdock.log'), 'w') as f:
                f.write(proc.stdout)

            out = root + '.sd' if os.path.isfile(root + '.sd') else root  # rbdock appends .sd to the -o root
            text = ''
            if os.path.isfile(out):
                with open(out) as f:
                    text = f.read()
            poses = {}
            for record in split_sd(text):
                poses.setdefault(sd_title(record), []).append(record)
            errors = RBDOCK_ERROR.findall(proc.stdout)  # [(error name, message line), ...]
            rejected = [msg or name for name, msg in errors if name == 'RBT_LIGAND_ERROR']
            if rejected:  # e.g. a ligand without the features of the mandatory pharmacophore restraints
                logger.warning('rbdock rejected %s: %s', sent[0] if len(sent) == 1 else f'{len(rejected)} ligand(s)',
                               '; '.join(rejected))
            if not poses:
                fatal = [f'{name}: {msg}' if msg else name for name, msg in errors if name != 'RBT_LIGAND_ERROR']
                tail = proc.stdout[-2000:]
                if fatal:  # receptor, cavity or protocol problem: rbdock exits 0, the request must fail instead
                    return {'status': 'error', 'error': f'rbdock failed ({fatal[0]}); check the receptor '
                                                        f'preparation. Output tail:\n{tail}'}
                if not rejected:
                    logger.warning('rbdock wrote no poses (exit %d). Output tail:\n%s', proc.returncode, tail)
                return {'status': 'ok', 'results': results}

            for name in payload:
                scored = []
                for record in poses.get(name, []):
                    try:
                        scored.append((float(sd_field(record, self.sort_field)),
                                       float(sd_field(record, self.score_field)), record))
                    except (TypeError, ValueError):
                        logger.warning('%s: a pose without numeric %s/%s fields was skipped',
                                       name, self.sort_field, self.score_field)
                if not scored:
                    logger.warning('%s: no poses from rbdock', name)
                    continue
                scored.sort(key=lambda t: t[0])  # stable: ties keep rbdock's order
                top = scored[:self.n_poses]
                results[name] = {
                    'docking_score': top[0][1],
                    'mol_block': pose_mol_block(top[0][2], name),
                    'raw_block': ''.join(t[2] for t in top),
                }
            return {'status': 'ok', 'results': results}
        except Exception as e:
            logger.exception('dock failed')
            return {'status': 'error', 'error': f'{type(e).__name__}: {e}'}
        finally:
            if self.keep_tmpdir:
                logger.info('docking files kept in %s', workdir)
            else:
                shutil.rmtree(workdir, ignore_errors=True)


def _emit(req_id, payload):
    sys.stdout.write(json.dumps({'id': req_id, 'payload': payload}) + '\n')
    sys.stdout.flush()


def run_server():
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    server = _RDockServer()
    # EasyDock stops the server with SIGTERM: turn it into a normal exit so the receptor folder is removed
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(0))

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue

        req_id = None
        try:
            request = json.loads(line)
            req_id = request.get('id')
            command = request.get('command')
            payload = request.get('payload') or {}

            if command == 'info':
                result = INFO.copy()
            elif command == 'init':
                result = server.init(payload)
            elif command == 'dock':
                result = server.dock(payload)
            else:
                result = {'status': 'error', 'error': f'Unknown command: {command}'}

        except Exception as e:
            result = {'status': 'error', 'error': str(e)}

        _emit(req_id, result)

    return 0


if __name__ == '__main__':
    raise SystemExit(run_server())
