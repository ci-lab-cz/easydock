import gzip
import logging
import os
import sqlite3
import sys
from typing import Optional, Dict, Union, Tuple, List

from rdkit import Chem

from easydock.auxiliary import count_input_structures
from easydock.database import get_protonation_arg_value, get_variables, _get_input_fname_from_setup


def get_pipeline_statistics(db_fname: str) -> Dict[str, Union[bool, int]]:
    """
    Return cumulative preparation and docking statistics from DB.
    """
    with sqlite3.connect(db_fname, timeout=90) as conn:
        cur = conn.cursor()

        protonation_enabled = get_protonation_arg_value(conn)

        try:
            input_structures_total = get_variables(conn, 'run_dock', ['input_structures_total'])['input_structures_total']
        except (sqlite3.OperationalError, KeyError) as e:
            logging.warning(f"Failed to read 'run_dock.input_structures_total' from variables table: {e}")
            input_structures_total = count_input_structures(_get_input_fname_from_setup(conn))

        successfully_read_structures = cur.execute(
            "SELECT COUNT(rowid) FROM mols WHERE stereo_id = 0"
        ).fetchone()[0]

        if input_structures_total is None:
            input_structures_total = successfully_read_structures

        total_generated_stereoisomers = cur.execute(
            "SELECT COUNT(rowid) FROM mols WHERE smi IS NOT NULL OR source_mol_block IS NOT NULL"
        ).fetchone()[0]

        if protonation_enabled:
            protonated_stereoisomers = cur.execute(
                "SELECT COUNT(rowid) FROM mols WHERE smi_protonated IS NOT NULL"
            ).fetchone()[0]
        else:
            protonated_stereoisomers = 0

        docked_stereoisomers = cur.execute(
            "SELECT COUNT(rowid) FROM mols WHERE docking_score IS NOT NULL"
        ).fetchone()[0]

    failed_stereoisomer_generation = get_stage_failures(db_fname, stage_name='stereoisomer_generation', output='count')
    failed_protonation = get_stage_failures(db_fname, stage_name='protonation', output='count')
    failed_docking = get_stage_failures(db_fname, stage_name='docking', output='count')

    return {
        'protonation_enabled': protonation_enabled,
        'input_structures': input_structures_total,
        'successfully_read_structures': successfully_read_structures,
        'generated_stereoisomers': total_generated_stereoisomers,
        'failed_stereoisomer_generation': failed_stereoisomer_generation,
        'protonated_stereoisomers': protonated_stereoisomers,
        'failed_protonation': failed_protonation,
        'docked_stereoisomers': docked_stereoisomers,
        'failed_docking': failed_docking
    }


def _get_error_log_fname(db_fname: str, stage_name: str) -> str:
    # stage_name = _normalize_stage_name(stage_name)
    base = os.path.abspath(db_fname)
    return f"{base}_{stage_name}_error.log"


def _read_smi_failures(input_fname: str, prefix: Optional[str] = None) -> Tuple[List[Tuple[str, str, int]], List[Tuple[str, str, int]]]:
    failed = []
    valid = []
    with open(input_fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            items = line.split('\t')
            smi = items[0].strip()
            if len(items) > 1 and items[1].strip():
                mol_id = items[1].strip()
            else:
                mol_id = smi
            if prefix:
                mol_id = f'{prefix}-{mol_id}'
            if Chem.MolFromSmiles(smi, sanitize=True) is None:
                failed.append((smi, mol_id, 0))
            else:
                valid.append((smi, mol_id, 0))
    return failed, valid


def _iter_sdf_blocks(input_fname: str):
    if input_fname.lower().endswith('.gz'):
        opener = lambda p: gzip.open(p, 'rt', encoding='utf-8', errors='ignore')
    else:
        opener = lambda p: open(p, 'rt', encoding='utf-8', errors='ignore')

    with opener(input_fname) as f:
        block = []
        for line in f:
            block.append(line)
            if line.strip() == '$$$$':
                yield ''.join(block)
                block = []
        if block:
            yield ''.join(block)


def _read_sdf_failures(input_fname: str, prefix: Optional[str] = None) -> Tuple[List[Tuple[str, str, int]], List[Tuple[str, str, int]]]:
    failed = []
    valid = []
    for i, block in enumerate(_iter_sdf_blocks(input_fname), 1):
        title = block.split('\n', 1)[0].strip()
        mol = Chem.MolFromMolBlock(block, sanitize=True, removeHs=False)
        if mol is None:
            mol_id = title if title else f'entry_{i}'
            if prefix:
                mol_id = f'{prefix}-{mol_id}'
            failed.append(('', mol_id, 0))
        else:
            if title:
                mol_id = title
            else:
                try:
                    mol_id = Chem.MolToSmiles(mol, isomericSmiles=True)
                except Exception:
                    mol_id = f'entry_{i}'
            if prefix:
                mol_id = f'{prefix}-{mol_id}'
            try:
                smi = Chem.MolToSmiles(mol, isomericSmiles=True)
            except Exception:
                smi = ''
            valid.append((smi, mol_id, 0))
    return failed, valid


def _get_source_parse_failures(db_fname: str, input_fname: Optional[str] = None, prefix: Optional[str] = None) -> List[Tuple[str, str, int]]:
    with sqlite3.connect(db_fname, timeout=90) as conn:
        if input_fname is None:
            input_fname = _get_input_fname_from_setup(conn)

        if not input_fname or not os.path.isfile(input_fname):
            logging.warning(f"Cannot build source parsing/storage error log: input file '{input_fname}' is unavailable.")
            return []

        lower_fname = input_fname.lower()
        if lower_fname.endswith('.smi') or lower_fname.endswith('.smiles'):
            failed, valid = _read_smi_failures(input_fname, prefix=prefix)
        elif lower_fname.endswith('.sdf') or lower_fname.endswith('.sdf.gz'):
            failed, valid = _read_sdf_failures(input_fname, prefix=prefix)
        else:
            logging.warning(f"Source parsing/storage error log is not implemented for input format: {input_fname}")
            return []

        stored_ids = {
            str(row[0]) for row in conn.execute("SELECT id FROM mols WHERE stereo_id = 0")
        }

        for smi, mol_id, stereo_id in valid:
            if str(mol_id) not in stored_ids:
                failed.append((smi, mol_id, stereo_id))

    # keep original order and remove duplicates
    seen = set()
    uniq = []
    for row in failed:
        if row not in seen:
            seen.add(row)
            uniq.append(row)
    return uniq


def _get_stage_failure_sql(stage_name: str, protonation_enabled: bool) -> Optional[str]:
    if stage_name == 'stereoisomer_generation':
        return """
            SELECT COALESCE(smi_input, smi, '') AS smi, id, stereo_id
            FROM mols
            WHERE stereo_id = 0
              AND smi IS NULL
              AND source_mol_block IS NULL
            ORDER BY id
        """
    if stage_name == 'protonation':
        if not protonation_enabled:
            return None
        return """
            SELECT COALESCE(smi, smi_input, '') AS smi, id, stereo_id
            FROM mols
            WHERE smi IS NOT NULL
              AND smi_protonated IS NULL
            ORDER BY id, stereo_id
        """
    if stage_name == 'docking':
        if protonation_enabled:
            return """
                SELECT COALESCE(smi, smi_input, smi_protonated, '') AS smi, id, stereo_id
                FROM mols
                WHERE (smi_protonated IS NOT NULL OR source_mol_block_protonated IS NOT NULL)
                  AND docking_score IS NULL
                ORDER BY id, stereo_id
            """
        return """
            SELECT COALESCE(smi, smi_input, '') AS smi, id, stereo_id
            FROM mols
            WHERE (smi IS NOT NULL OR source_mol_block IS NOT NULL)
              AND docking_score IS NULL
            ORDER BY id, stereo_id
        """
    return None


def get_stage_failures(db_fname: str,
                       stage_name: str,
                       input_fname: Optional[str] = None,
                       prefix: Optional[str] = None,
                       output: str = 'records') -> Union[List[Tuple[str, str, int]], int]:
    """
    Return failed records for a given stage.
    :param db_fname:
    :param stage_name:
    :param input_fname:
    :param prefix:
    :param output: 'records' for [(smi, id, stereo_id)], 'count' for integer count
    """
    # stage_name = _normalize_stage_name(stage_name)
    if output not in ('records', 'count'):
        raise ValueError(f"Unknown output mode: {output}")

    if stage_name == 'input_parsing':
        records = _get_source_parse_failures(db_fname, input_fname=input_fname, prefix=prefix)
        if output == 'count':
            return len(records)
        return records

    with sqlite3.connect(db_fname, timeout=90) as conn:
        cur = conn.cursor()
        protonation_enabled = get_protonation_arg_value(conn)
        sql = _get_stage_failure_sql(stage_name, protonation_enabled)
        if sql is None:
            return 0 if output == 'count' else []

        if output == 'count':
            count_sql = f"SELECT COUNT(*) FROM ({sql}) AS q"
            return int(cur.execute(count_sql).fetchone()[0])

        rows = cur.execute(sql).fetchall()

    return [(smi if smi is not None else '', mol_id, stereo_id) for smi, mol_id, stereo_id in rows]


def write_stage_error_log(db_fname: str, stage_name: str, input_fname: Optional[str] = None, prefix: Optional[str] = None) -> Tuple[str, int]:
    """
    Create a stage error log file with columns: smi, id, stereo_id
    """
    records = get_stage_failures(db_fname, stage_name=stage_name, input_fname=input_fname, prefix=prefix, output='records')

    if records:
        fname = _get_error_log_fname(db_fname, stage_name)

        with open(fname, 'wt', encoding='utf-8') as f:
            f.write('smi\tid\tstereo_id\n')
            for smi, mol_id, stereo_id in records:
                f.write(f'{smi}\t{mol_id}\t{stereo_id}\n')
    else:
        fname = None

    return fname, len(records)


def report_error_log_file(stage_name: str, fname: str, nfailed: int):
    if fname and nfailed > 0:
        message = f'{stage_name} error log: {fname} ({nfailed} records)'
        logging.info(message)
        sys.stderr.write(message + '\n')
        sys.stderr.flush()
