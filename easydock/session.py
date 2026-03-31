"""
Session save/restore utilities for EasyDock.

Stores command-line args and config file contents (including all referenced text files)
in the database `setup` table, and reconstructs them on restart.

Setup table schema:
    CREATE TABLE IF NOT EXISTS setup (key TEXT PRIMARY KEY, content TEXT)

Key conventions:
    'args'             — JSON dump of argparse namespace
    'config'           — raw text of config.yml
    'file:<path>'      — content of a text file at dotted key path in the config YAML
                         e.g. 'file:protein', 'file:init_server.protein'
"""

import json
import logging
import os
import sqlite3
import tempfile

import yaml

from easydock.auxiliary import expand_path

logger = logging.getLogger(__name__)

DB_USER_VERSION = 1


# ── internal helpers ──────────────────────────────────────────────────────────

def _upsert(conn, key, content):
    conn.execute(
        'INSERT INTO setup (key, content) VALUES (?, ?) '
        'ON CONFLICT(key) DO UPDATE SET content = excluded.content',
        (key, content)
    )
    conn.commit()


def _is_text_file(path):
    """Return True if the file appears to be plain text (UTF-8, no null bytes)."""
    try:
        with open(path, 'rb') as f:
            chunk = f.read(8192)
        if b'\x00' in chunk:
            return False
        chunk.decode('utf-8')
        return True
    except (UnicodeDecodeError, IOError):
        return False


def _extract_file_paths(obj, prefix=''):
    """
    Recursively walk a nested dict/list and return [(dotted_key_path, abs_file_path), ...]
    for all string values that point to existing, readable text files.
    """
    results = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            child_key = f'{prefix}.{k}' if prefix else str(k)
            results.extend(_extract_file_paths(v, child_key))
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            results.extend(_extract_file_paths(v, f'{prefix}.{i}'))
    elif isinstance(obj, str):
        p = os.path.abspath(expand_path(obj))
        if os.path.isfile(p) and _is_text_file(p):
            results.append((prefix, p))
    return results


def _get_nested(d, dotted_key):
    """Return d[k1][k2]...[kn] for a dotted key, or None if any key is missing."""
    parts = dotted_key.split('.')
    for k in parts:
        if isinstance(d, dict) and k in d:
            d = d[k]
        elif isinstance(d, (list, tuple)):
            try:
                d = d[int(k)]
            except (ValueError, IndexError):
                return None
        else:
            return None
    return d


def _set_nested(d, dotted_key, value):
    """Set d[k1][k2]...[kn] = value for a dotted key. Silently skips missing paths."""
    parts = dotted_key.split('.')
    for k in parts[:-1]:
        if isinstance(d, dict) and k in d:
            d = d[k]
        elif isinstance(d, (list, tuple)):
            try:
                d = d[int(k)]
            except (ValueError, IndexError):
                return
        else:
            return
    last = parts[-1]
    if isinstance(d, dict):
        d[last] = value
    elif isinstance(d, (list, tuple)):
        try:
            d[int(last)] = value
        except (ValueError, IndexError):
            pass


# ── public API ────────────────────────────────────────────────────────────────

def create_setup_table(conn):
    """Create the setup key/value table if it does not exist."""
    conn.execute("""
        CREATE TABLE IF NOT EXISTS setup (
            key     TEXT PRIMARY KEY,
            content TEXT
        )
    """)
    conn.commit()


def save_session_args(conn, args_dict):
    """Store argparse namespace dict as JSON under the 'args' key."""
    _upsert(conn, 'args', json.dumps(args_dict, indent=2))


def save_session_config(conn, config_path):
    """
    Store the config.yml text and the content of every text file referenced inside it.

    :param conn: open SQLite connection
    :param config_path: path to config YAML file
    """
    config_path = os.path.abspath(expand_path(config_path))
    with open(config_path) as f:
        raw = f.read()
    _upsert(conn, 'config', raw)

    config_dict = yaml.safe_load(raw) or {}
    if isinstance(config_dict, dict):
        for dotted_path, abs_path in _extract_file_paths(config_dict):
            with open(abs_path) as f:
                file_content = f.read()
            _upsert(conn, f'file:{dotted_path}', file_content)


def restore_session(conn, allowed_override_args=(), supplied_args=(), tmpdir=None):
    """
    Restore args dict and recreate temp files from the setup table.

    :param conn: open SQLite connection
    :param allowed_override_args: set/list of arg names that CLI values may override
    :param supplied_args: set/list of arg names actually typed on the CLI
                          (from get_supplied_args); only these can override stored values
    :param tmpdir: directory for temp files; uses system default if None
    :return: (args_dict, tmpfiles) where tmpfiles is a list of paths to delete on exit
    """
    ver = conn.execute('PRAGMA user_version').fetchone()[0]
    if ver < DB_USER_VERSION:
        raise RuntimeError(
            f'Database user_version={ver} is older than the required {DB_USER_VERSION}. '
            'Run easydock_migrate to upgrade it.'
        )

    rows = dict(conn.execute('SELECT key, content FROM setup').fetchall())

    if 'args' not in rows:
        raise RuntimeError("'args' key not found in setup table — database may be corrupted.")

    # 1. Restore args
    args_dict = json.loads(rows['args'])

    # 2. Apply CLI overrides: remove stored values for args that were explicitly
    #    supplied on the command line AND are in the allowed-override list.
    #    supplied_args comes from get_supplied_args() which scans sys.argv directly,
    #    so only flags the user actually typed are included — never argparse defaults.
    allowed = set(allowed_override_args)
    supplied = set(supplied_args)
    for arg in supplied & allowed:
        args_dict.pop(arg, None)

    tmpfiles = []
    backup_tempdir = tempfile.tempdir
    if tmpdir:
        tempfile.tempdir = tmpdir
    try:
        # 3. Rebuild config YAML with updated (temp) file paths
        if 'config' in rows:
            config_dict = yaml.safe_load(rows['config']) or {}
            file_entries = [(k[len('file:'):], v) for k, v in rows.items()
                            if k.startswith('file:')]
            for dotted_path, content in file_entries:
                orig = _get_nested(config_dict, dotted_path)
                if orig is None:
                    logger.warning(
                        'session restore: config has no key %r — skipping temp file.', dotted_path
                    )
                    continue
                suffix = os.path.splitext(orig)[1] or '.tmp'
                fd, tmp_path = tempfile.mkstemp(suffix=suffix, text=True)
                os.close(fd)
                with open(tmp_path, 'wt') as fh:
                    fh.write(content)
                tmpfiles.append(tmp_path)
                _set_nested(config_dict, dotted_path, tmp_path)

            fd, cfg_tmp = tempfile.mkstemp(suffix='.yml', text=True)
            os.close(fd)
            with open(cfg_tmp, 'wt') as fh:
                yaml.safe_dump(config_dict, fh)
            tmpfiles.append(cfg_tmp)
            args_dict['config'] = cfg_tmp

    except Exception:
        for f in tmpfiles:
            try:
                os.unlink(f)
            except OSError:
                pass
        raise
    finally:
        tempfile.tempdir = backup_tempdir

    return args_dict, tmpfiles


def migrate_setup_table(db_fname):
    """
    Migrate a database with the old dynamic-column setup table to the new key/value schema.

    Old schema:  setup(yaml TEXT, config TEXT, protein TEXT, protein_setup TEXT, ...)
    New schema:  setup(key TEXT PRIMARY KEY, content TEXT)
    """
    with sqlite3.connect(db_fname, timeout=90) as conn:
        ver = conn.execute('PRAGMA user_version').fetchone()[0]
        if ver >= DB_USER_VERSION:
            logger.info('migrate_setup_table: database already at version %d, skipping.', ver)
            return
        if ver != 0:
            raise RuntimeError(
                f'migrate_setup_table: expected version 0, got {ver}. '
                'This migration only handles version 0 → 1.'
            )

        # Read the single row
        res = conn.execute('SELECT * FROM setup').fetchone()
        if res is None:
            raise RuntimeError('setup table is empty; nothing to migrate.')
        col_names = [d[0] for d in conn.execute('SELECT * FROM setup').description]
        old_row = dict(zip(col_names, res))

        # Drop and recreate
        conn.execute('DROP TABLE setup')
        create_setup_table(conn)

        # Migrate 'yaml' → 'args' (convert from YAML to JSON)
        if old_row.get('yaml'):
            args_dict = yaml.safe_load(old_row['yaml']) or {}
            _upsert(conn, 'args', json.dumps(args_dict, indent=2))

        # Migrate 'config' → 'config'
        if old_row.get('config'):
            _upsert(conn, 'config', old_row['config'])

        # Migrate remaining columns (file contents) → 'file:<colname>'
        for col in col_names:
            if col in ('yaml', 'config') or old_row.get(col) is None:
                continue
            _upsert(conn, f'file:{col}', old_row[col])

        # Set new user_version
        conn.execute(f'PRAGMA user_version = {DB_USER_VERSION}')
        conn.commit()

    logger.info('migrate_setup_table: migration complete for %s', db_fname)
