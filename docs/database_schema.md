# Database Schema

This page describes the internal SQLite database structure used by EasyDock. It is intended for developers extending the pipeline or reading results programmatically.

---

## Overview

EasyDock stores all molecule data, docking results, and session configuration in a single **SQLite** file (the `-o` / `--output` argument). The database is opened in **WAL (Write-Ahead Logging)** mode to allow concurrent reads from result-retrieval tools while docking is in progress.

The schema consists of six tables:

| Table | Purpose |
|---|---|
| [`mols`](#mols) | Molecule structures and docking results — the main data table |
| [`setup`](#setup) | Session configuration for resuming interrupted runs |
| [`variables`](#variables) | Module-level metadata and run-time flags |
| [`plif_names`](#plif_names) | Lookup table for protein-ligand interaction contact types |
| [`plif_res`](#plif_res) | Per-molecule, per-pose interaction fingerprint results |
| [`bust`](#bust) | PoseBusters validation results per pose |

---

## `mols`

The central table. One row per enumerated stereoisomer of each input molecule. The primary key is `(id, stereo_id)`.

| Column | Type | Default | Description                                                                                                                         |
|---|---|---|-------------------------------------------------------------------------------------------------------------------------------------|
| `id` | TEXT | — | Molecule identifier — taken from the input<br>file title or derived from SMILES.                                                    |
| `stereo_id` | INTEGER | `0` | Stereoisomer index. `0` for the first (or only)<br>isomer; incremented for each additional one.                                     |
| `smi_input` | TEXT | NULL | Raw SMILES exactly as read from the input file.<br>Only set for SMILES inputs; NULL for SDF/PDBQT.                                  |
| `smi` | TEXT | NULL | Canonical SMILES after salt stripping andvz<br>stereoisomer assignment. Used as the docking<br>input for 2D structures.             |
| `smi_protonated` | TEXT | NULL | Protonated canonical SMILES at the target pH.<br>Set only when `--protonation` is used.                                             |
| `source_mol_block_input` | TEXT | NULL | Original mol block from the input SDF/PDBQT file.<br>Set only for 3D structure inputs.                                              |
| `source_mol_block` | TEXT | NULL | Cleaned mol block (salts stripped, 2D/3D structure)<br>prepared for docking. Set for 3D inputs.                                     |
| `source_mol_block_protonated` | TEXT | NULL | Mol block after protonation, inheriting 3D<br>coordinates from `source_mol_block`.<br>Set when protonation is applied to<br>a 3D input. |
| `docking_score` | REAL | NULL | Best docking score in kcal/mol. NULL until<br>docking completes for this molecule.                                                  |
| `raw_block` | TEXT | NULL | All docking poses in their native format<br>(PDBQT or SDF, controlled by the `raw_format`<br>variable). NULL until docking completes.  |
| `mol_block` | TEXT | NULL | Best docking pose as an MDL mol block.<br>NULL until docking completes.                                                             |
| `dock_time` | REAL | NULL | Wall-clock docking time for this molecule in<br>seconds.                                                                               |
| `time` | TEXT | NULL | ISO timestamp of the last write to this row<br>(`datetime(current_timestamp, 'localtime')`).<br>Updated on every `UPDATE`.          |

### Primary Key

```sql
PRIMARY KEY (id, stereo_id)
```

### Processing Stages

Each stage selects its work by checking which columns are still NULL:

| Stage | Selects rows where | Writes columns |
|---|---|---|
| **Input parsing** | (new rows) | `id`, `stereo_id`, `smi_input`, `smi`,<br>`source_mol_block_input`,<br>`source_mol_block` |
| **Protonation** | `smi IS NOT NULL`<br>`AND smi_protonated IS NULL`<br>`AND docking_score IS NULL` | `smi_protonated`,<br>`source_mol_block_protonated` |
| **Docking** (with protonation) | `docking_score IS NULL`<br>`AND (smi_protonated IS NOT NULL`<br>`OR source_mol_block_protonated IS NOT NULL)` | `docking_score`, `raw_block`,<br>`mol_block`, `dock_time` |
| **Docking** (without protonation) | `docking_score IS NULL`<br>`AND (smi IS NOT NULL`<br>`OR source_mol_block IS NOT NULL)` | `docking_score`, `raw_block`,<br>`mol_block`, `dock_time` |

### Notes

- Rows with `smi IS NULL AND source_mol_block IS NULL AND stereo_id = 0` represent molecules that failed stereoisomer generation and are silently skipped by all downstream stages.
- `mol_block` stores only the best pose; all poses are in `raw_block`. Use `raw_block` when re-ranking or re-scoring poses.
- The `raw_block` format (PDBQT vs SDF) is determined by the `raw_format` variable in the [`variables`](#variables) table.

---

## `setup`

Stores the full session configuration so interrupted runs can be resumed with `--input` and `--output` only. One row per key.

| Column | Type | Description |
|---|---|---|
| `key` | TEXT (PRIMARY KEY) | Configuration key (see below). |
| `content` | TEXT | Configuration value — raw text, YAML, or JSON depending on the key. |

### Keys

| Key | Content |
|---|---|
| `args` | JSON-serialised `argparse` namespace — all CLI arguments as passed to the original run. |
| `config` | Raw YAML text of the docking config file. |
| `file:<path>` | Content of a text file referenced inside the config YAML. `<path>` is the dotted YAML key<br>(e.g., `file:protein`, `file:init_server.protein`).<br>Stored so the session can be fully reconstructed without the original files. |

---

## `variables`

A generic key-value store for module-level metadata. Allows any module to persist flags or settings without a schema change. Composite primary key `(module, name)`.

| Column | Type | Description |
|---|---|---|
| `module` | TEXT (NOT NULL) | Module that owns this variable (e.g., `run_dock`, `database`, `easydock_bust`). |
| `name` | TEXT (NOT NULL) | Variable name. |
| `value_type` | TEXT | Type hint for deserialisation: `'int'`, `'float'`, `'str'`, or `'json'`. |
| `value` | TEXT | Stringified value. |
| `updated_at` | TIMESTAMP | Automatically set to `CURRENT_TIMESTAMP` on insert/update. |

### Well-Known Variables

| Module | Name | Type | Values | Meaning |
|---|---|---|---|---|
| `database` | `raw_format` | str | `'pdbqt'`, `'sdf'` | Format of pose data stored in `mols.raw_block`.<br>Defaults to `'pdbqt'`; set to `'sdf'` for<br>server-based docking programs that return SDF. |
| `run_dock` | `input_structures_total` | int | positive integer | Total number of input structures counted at startup —<br>used for progress reporting. |
| `easydock_bust` | `bust_protein` | str | PDB text | Protein PDB block with explicit hydrogens added,<br>used by PoseBusters. |
| `easydock_plif` | `plif_protein` | str | PDB text | Protein PDB block with explicit hydrogens added,<br>used by ProLIF. |

---

## `plif_names`

Lookup table for protein-ligand interaction contact types, populated by `easydock_plif`. Avoids repeating long contact-name strings in `plif_res`.

| Column | Type | Description |
|---|---|---|
| `plif_id` | INTEGER (PRIMARY KEY AUTOINCREMENT) | Surrogate key for the contact type. |
| `contact_name` | TEXT (UNIQUE) | Human-readable contact identifier,<br>e.g., `'GLU80.A.HBDonor'`.<br>Format: `<residue>.<chain>.<interaction_type>`. |

---

## `plif_res`

One row per detected contact, per pose, per molecule. References `mols` by `rowid` and `plif_names` by `plif_id`.

| Column | Type | Description |
|---|---|---|
| `mols_rowid` | INTEGER | `rowid` of the corresponding row in `mols`. |
| `pose` | INTEGER | 1-based pose index (`1` = best docking pose). |
| `plif_id` | INTEGER | Foreign key into `plif_names`. |

**Unique constraint:** `(mols_rowid, pose, plif_id)`

To retrieve contacts for a molecule:

```sql
SELECT pn.contact_name
FROM plif_res pr
JOIN plif_names pn ON pr.plif_id = pn.plif_id
JOIN mols m ON pr.mols_rowid = m.rowid
WHERE m.id = 'CHEMBL12345' AND m.stereo_id = 0 AND pr.pose = 1;
```

---

## `bust`

One row per tested pose per molecule. Populated by `easydock_bust`.

| Column | Type | Description |
|---|---|---|
| `mols_rowid` | INTEGER | `rowid` of the corresponding row in `mols`. |
| `pose` | INTEGER | 1-based pose index. |
| `result` | INTEGER | `1` if the pose passes all PoseBusters checks; `0` otherwise. |

**Unique constraint:** `(mols_rowid, pose)`

---

## Schema Version

EasyDock tracks schema version via SQLite's built-in `PRAGMA user_version`. The current version is **1**.

On database open, `session.py` checks the version and runs any necessary migrations:

| Migration | From → To | What changes |
|---|---|---|
| `migrate_setup_table()` | 0 → 1 | Converts the old `setup` table (dynamic columns, one row)<br>to the current key-value schema. |
| `migrate_pdb_block_to_raw_block()` | legacy → current | Renames the `pdb_block` column in `mols` to `raw_block`<br>if the old name is still present.<br>Run automatically when `raw_format` variable is absent. |

The user is prompted before any migration is applied (120-second timeout; auto-abort on no response).

---

## Concurrency Model

- **WAL mode** is enabled on every open connection: `PRAGMA journal_mode=WAL`.
- During docking, the main process writes results via `update_db()` while worker threads read molecules via `MolQueue`.
- `MolQueue` opens a **separate read connection** per batch, closes it between fetches, and tracks position by `rowid` — allowing the writer to insert/update rows concurrently without blocking readers.
- Result-retrieval tools (`get_sdf`, `easydock_bust`, `easydock_plif`) may safely open the database read-only while a run is in progress.
