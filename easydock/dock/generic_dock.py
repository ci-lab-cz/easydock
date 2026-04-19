import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import timeit

import yaml
from rdkit import Chem

from easydock.auxiliary import expand_path, resolve_path
from easydock.dock.preparation_for_docking import ligand_preparation, pdbqt2molblock

logger = logging.getLogger(__name__)

_IN_FORMAT_SUFFIX = {
    'pdbqt': '.pdbqt',
    'smiles': '.smi',
    'mol': '.sdf',
}

_OUT_FORMAT_SUFFIX = {
    'sdf': '.sdf',
    'pdbqt': '.pdbqt',
    'pdb': '.pdb',
}


def __parse_config(config_fname):
    with open(config_fname) as f:
        config = yaml.safe_load(f) or {}
    config_dir = os.path.dirname(os.path.abspath(config_fname))
    if 'script_file' in config:
        config['script_file'] = resolve_path(str(config['script_file']), config_dir)
    if 'env' in config:
        config['env'] = resolve_path(str(config['env']), config_dir)
    program_args = config.get('program_args') or {}
    for key, value in program_args.items():
        if isinstance(value, str):
            program_args[key] = resolve_path(value, config_dir)
    return config


def _find_conda_executable():
    """Return the first available conda-family executable (conda, mamba, micromamba)."""
    for name in ('conda', 'mamba', 'micromamba'):
        path = shutil.which(name)
        if path:
            return path
    raise RuntimeError(
        "No conda/mamba executable found in PATH. "
        "Install conda or mamba, or specify 'env' as the full path to the environment directory."
    )


def _python_prefix(script_file, env):
    """Return the command prefix [interpreter, script_file] for a Python script.

    env may be:
      - None              → use the current Python interpreter (sys.executable)
      - a directory path  → use {env}/bin/python directly (venv or conda env by path)
      - an env name       → use ``conda run -n {env}`` (auto-detects conda/mamba/micromamba)
    """
    if env is None:
        return [sys.executable, script_file]

    env_str = str(env)
    env_expanded = expand_path(env_str)

    if os.path.isdir(env_expanded):
        # Explicit path to a venv or conda environment directory
        python = os.path.join(env_expanded, 'bin', 'python')
        if not os.path.isfile(python):
            python = os.path.join(env_expanded, 'Scripts', 'python.exe')  # Windows
        if not os.path.isfile(python):
            raise RuntimeError(f"No Python interpreter found in environment directory: {env_expanded}")
        return [python, script_file]

    # Treat as a conda environment name
    conda = _find_conda_executable()
    # --no-capture-output ensures the script's stdout/stderr flow through to our subprocess
    return [conda, 'run', '--no-capture-output', '-n', env_str, 'python', script_file]


def _build_base_cmd(config):
    script_file = config['script_file']
    env = config.get('env')

    if script_file.endswith('.py'):
        cmd = _python_prefix(script_file, env)
    else:
        if env is not None:
            logger.warning(
                "'env' is set but script_file %r does not end in .py; 'env' will be ignored",
                script_file,
            )
        cmd = shlex.split(script_file)

    program_args = config.get('program_args') or {}
    for key, value in program_args.items():
        if key == 'extra_args':
            continue
        if value is not None:
            cmd += [f'--{key}', str(value)]
    extra_args = program_args.get('extra_args', '')
    if extra_args:
        cmd += shlex.split(str(extra_args))
    return cmd


def _prepare_ligand(mol, ligand_in_format, ring_sample):
    if ligand_in_format == 'pdbqt':
        return ligand_preparation(mol, boron_replacement=True, ring_sample=ring_sample)
    if ligand_in_format == 'smiles':
        return [Chem.MolToSmiles(mol, isomericSmiles=True)]
    if ligand_in_format == 'mol':
        return [Chem.MolToMolBlock(mol)]
    raise ValueError(f'Unsupported ligand_in_format: {ligand_in_format!r}')


def _parse_score(raw_block, out_format, parse_score_cfg):
    """Extract score from the first pose in raw_block (poses are pre-sorted by program)."""
    if not raw_block or not raw_block.strip():
        return None
    score_field = parse_score_cfg.get('score_field')
    score_regex = parse_score_cfg.get('score_regex')

    if out_format == 'sdf':
        first = raw_block.split('$$$$')[0]
        if score_field:
            m = re.search(
                rf'>\s*<{re.escape(score_field)}>\s*\n\s*(-?[\d.]+(?:[eE][+-]?\d+)?)',
                first)
            if m:
                try:
                    return float(m.group(1))
                except ValueError:
                    pass

    elif out_format in ('pdbqt', 'pdb'):
        parts = re.split(r'^MODEL\b', raw_block, maxsplit=1, flags=re.MULTILINE)
        first = parts[1] if len(parts) > 1 else raw_block
        if score_regex:
            m = re.search(score_regex, first)
            if m:
                try:
                    return float(m.group(1))
                except (ValueError, IndexError):
                    pass

    return None


def _build_mol_block(raw_block, out_format, mol, mol_id):
    if out_format == 'sdf':
        return raw_block.split('$$$$')[0]
    if out_format == 'pdbqt':
        parts = re.split(r'^MODEL\b', raw_block, maxsplit=1, flags=re.MULTILINE)
        first_model = parts[1] if len(parts) > 1 else raw_block
        try:
            return pdbqt2molblock(first_model, mol, mol_id)
        except Exception:
            logger.exception('Failed to convert pdbqt pose to mol_block for %s', mol_id)
            return None
    if out_format == 'pdb':
        parts = re.split(r'^MODEL\b', raw_block, maxsplit=1, flags=re.MULTILINE)
        return (parts[1] if len(parts) > 1 else raw_block).strip()
    return None


def _run_one(ligand_str, base_cmd, config):
    """Run the binary for one ligand. Returns (score, raw_block) or (None, None)."""
    out_format = config['ligand_out_format']
    in_format = config['ligand_in_format']
    parse_score_cfg = config.get('parse_score') or {}
    input_arg_name = config.get('input_arg_name')
    output_arg_name = config.get('output_arg_name')

    use_stdin = input_arg_name is None
    use_stdout = output_arg_name is None

    in_suffix = _IN_FORMAT_SUFFIX.get(in_format, '.mol')
    out_suffix = _OUT_FORMAT_SUFFIX.get(out_format, '.out')

    in_fd = in_fname = out_fd = out_fname = None
    try:
        if not use_stdin:
            in_fd, in_fname = tempfile.mkstemp(suffix=in_suffix)
            os.close(in_fd)
            in_fd = None
            with open(in_fname, 'w') as f:
                f.write(ligand_str)

        if not use_stdout:
            out_fd, out_fname = tempfile.mkstemp(suffix=out_suffix)
            os.close(out_fd)
            out_fd = None

        cmd = list(base_cmd)
        if not use_stdin:
            cmd += [input_arg_name, in_fname]
        if not use_stdout:
            cmd += [output_arg_name, out_fname]

        logger.debug('Running: %s', ' '.join(str(x) for x in cmd))
        proc = subprocess.run(
            cmd,
            input=ligand_str if use_stdin else None,
            capture_output=use_stdout,
            text=True,
        )

        if proc.stderr:
            logger.debug('stderr:\n%s', proc.stderr)
        if proc.returncode != 0:
            logger.warning(
                'Binary exited with code %d\nstdout: %s\nstderr: %s',
                proc.returncode,
                proc.stdout if use_stdout else '',
                proc.stderr,
            )
            return None, None

        if use_stdout:
            raw_block = proc.stdout
        else:
            with open(out_fname) as f:
                raw_block = f.read()

        score = _parse_score(raw_block, out_format, parse_score_cfg)
        return score, raw_block

    finally:
        if in_fd is not None:
            os.close(in_fd)
        if out_fd is not None:
            os.close(out_fd)
        if in_fname and os.path.exists(in_fname):
            os.unlink(in_fname)
        if out_fname and os.path.exists(out_fname):
            os.unlink(out_fname)


def mol_dock(mol, config, ring_sample=False):
    """
    Dock a single molecule using an arbitrary external binary.

    :param mol: RDKit Mol with _Name property
    :param config: path to YAML config file
    :param ring_sample: sample saturated ring conformers and dock each
    :return: (mol_id, result_dict) or (mol_id, None) on failure
    """
    config = __parse_config(config)
    mol_id = mol.GetProp('_Name')
    out_format = config['ligand_out_format']
    score_mode = config.get('score_mode', 'min')
    base_cmd = _build_base_cmd(config)

    ligand_list = _prepare_ligand(mol, config['ligand_in_format'], ring_sample)
    if not ligand_list:
        return mol_id, None

    results = []
    start_time = timeit.default_timer()

    for ligand_str in ligand_list:
        score, raw_block = _run_one(ligand_str, base_cmd, config)
        if score is None:
            logger.warning('No score extracted for %s', mol_id)
            continue
        mol_block = _build_mol_block(raw_block, out_format, mol, mol_id)
        if mol_block:
            results.append({'docking_score': score, 'raw_block': raw_block, 'mol_block': mol_block})

    dock_time = round(timeit.default_timer() - start_time, 1)

    if not results:
        return mol_id, None

    choose = min if score_mode == 'min' else max
    output = choose(results, key=lambda x: x['docking_score'])
    output['dock_time'] = dock_time
    return mol_id, output
