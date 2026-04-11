import json
import logging
import os
import re
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)

INFO = {
    'batch_size': 1,
    'ligand_in_format': 'pdbqt',
    'ligand_out_format': 'pdbqt',
    'score_mode': 'min',
}

PROGRAMS = {
    # GPU programs (OpenCL-accelerated)
    'vina-gpu': {
        'binary': '/opt/vinagpu/vina-gpu/AutoDock-Vina-GPU-2-1',
        'opencl_path': '/opt/vinagpu/vina-gpu',
        'gpu': True,
    },
    'qvina-gpu': {
        'binary': '/opt/vinagpu/qvina-gpu/QuickVina2-GPU-2-1',
        'opencl_path': '/opt/vinagpu/qvina-gpu',
        'gpu': True,
    },
    'qvinaw-gpu': {
        'binary': '/opt/vinagpu/qvinaw-gpu/QuickVina-W-GPU-2-1',
        'opencl_path': '/opt/vinagpu/qvinaw-gpu',
        'gpu': True,
    },
    # CPU programs from QVina/AutoDock-Vina repos
    'vina': {
        'binary': '/opt/vinagpu/vina/vina',
        'gpu': False,
    },
    'qvina': {
        'binary': '/opt/vinagpu/qvina/qvina2.1',
        'gpu': False,
    },
    'qvinaw': {
        'binary': '/opt/vinagpu/qvinaw/qvina-w',
        'gpu': False,
    },
}


def _expand_path(value):
    if value is None:
        return None
    return str(Path(value).expanduser().resolve())


def _parse_score(pdbqt_text):
    match = re.search(r'REMARK VINA RESULT:\s+(-?[\d.]+)', pdbqt_text)
    if match:
        return round(float(match.group(1)), 3)
    return None


class _VinaGPUServer:
    def __init__(self):
        self.initialized = False
        self.binary = None
        self.is_gpu = None
        self.protein = None
        self.protein_setup = None
        # GPU-only
        self.opencl_binary_path = None
        self.thread = 8000
        # CPU-only
        self.exhaustiveness = 8
        self.cpu = None
        # common
        self.n_poses = 9
        self.seed = None
        self.extra_args = []

    def init(self, payload):
        program = payload.get('program')
        if program not in PROGRAMS:
            return {'status': 'error',
                    'error': f"Unknown program {program!r}. Choose from: {list(PROGRAMS)}"}

        prog = PROGRAMS[program]
        self.binary = prog['binary']
        self.is_gpu = prog['gpu']
        self.protein = _expand_path(payload.get('protein'))
        self.protein_setup = _expand_path(payload.get('protein_setup'))
        self.n_poses = int(payload.get('n_poses', self.n_poses))
        self.seed = payload.get('seed')
        extra_args_raw = payload.get('extra_args', [])
        if isinstance(extra_args_raw, str):
            self.extra_args = shlex.split(extra_args_raw)
        else:
            self.extra_args = [str(a) for a in extra_args_raw]

        if self.is_gpu:
            self.opencl_binary_path = payload.get('opencl_binary_path') or prog['opencl_path']
            self.thread = int(payload.get('thread', self.thread))
        else:
            self.exhaustiveness = int(payload.get('exhaustiveness', self.exhaustiveness))
            cpu = payload.get('cpu')
            self.cpu = int(cpu) if cpu is not None else None

        if not self.protein:
            return {'status': 'error', 'error': "'protein' is required"}
        if not self.protein_setup:
            return {'status': 'error', 'error': "'protein_setup' is required"}

        self.initialized = True
        return {'status': 'ok'}

    def dock(self, payload):
        if not self.initialized:
            return {'status': 'error', 'error': 'Server not initialized'}

        if not isinstance(payload, dict):
            return {'status': 'error',
                    'error': f'payload must be a dict of {{mol_id: pdbqt_string}}. Got: {type(payload)}'}

        results = {}
        for mol_id, pdbqt_string in payload.items():
            if not pdbqt_string:
                logger.warning('Empty pdbqt for %s', mol_id)
                results[mol_id] = None
                continue

            out_fd, out_fname = tempfile.mkstemp(suffix='_output.pdbqt')
            lig_fd, lig_fname = tempfile.mkstemp(suffix='_ligand.pdbqt')
            try:
                with open(lig_fname, 'w') as f:
                    f.write(pdbqt_string)

                cmd = [
                    self.binary,
                    '--receptor', self.protein,
                    '--ligand', lig_fname,
                    '--out', out_fname,
                    '--config', self.protein_setup,
                    '--num_modes', str(self.n_poses),
                ]
                if self.is_gpu:
                    cmd += [
                        '--opencl_binary_path', self.opencl_binary_path,
                        '--thread', str(self.thread),
                    ]
                else:
                    cmd += ['--exhaustiveness', str(self.exhaustiveness)]
                    if self.cpu is not None:
                        cmd += ['--cpu', str(self.cpu)]
                if self.seed is not None:
                    cmd += ['--seed', str(self.seed)]
                cmd += self.extra_args

                logger.debug('Running: %s', ' '.join(cmd))
                proc = subprocess.run(cmd, capture_output=True, text=True)
                if proc.stdout:
                    logger.debug('stdout for %s:\n%s', mol_id, proc.stdout)
                if proc.stderr:
                    logger.info('stderr for %s:\n%s', mol_id, proc.stderr)
                if proc.returncode != 0:
                    logger.warning(
                        'Binary exited with code %d for %s\nstdout: %s\nstderr: %s',
                        proc.returncode, mol_id, proc.stdout, proc.stderr,
                    )
                    results[mol_id] = None
                    continue

                with open(out_fname) as f:
                    raw_block = f.read()

                score = _parse_score(raw_block)
                if score is None:
                    logger.warning(
                        'No score found in output for %s; output file content:\n%s',
                        mol_id, raw_block or '<empty>',
                    )
                    results[mol_id] = None
                else:
                    results[mol_id] = {
                        'docking_score': score,
                        'raw_block': raw_block,
                        'mol_block': None,
                    }

            except Exception:
                logger.exception('Unexpected error docking %s', mol_id)
                results[mol_id] = None
            finally:
                os.close(out_fd)
                os.close(lig_fd)
                if os.path.exists(lig_fname):
                    os.unlink(lig_fname)
                if os.path.exists(out_fname):
                    os.unlink(out_fname)

        return {'status': 'ok', 'results': results}


def _emit(req_id, payload):
    sys.stdout.write(json.dumps({'id': req_id, 'payload': payload}) + '\n')
    sys.stdout.flush()


def run_server():
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    server = _VinaGPUServer()

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
