# Generic Binary Docking

`--program binary` is a generic wrapper that can run **any external docking binary or Python
script** through a single YAML config file.  No code changes are needed — the config drives
everything: which executable to call, how to pass the ligand and receive the output, and how to
extract the score.

```bash
easydock -i input.smi -o output.db --program binary --config config.yml -c 4 --sdf
```

## Config reference

| Key | Required | Description |
|-----|----------|-------------|
| `script_file` | yes | Executable or command string. `~` and environment variables are expanded first, then the result is split on whitespace, so `apptainer run gnina.sif gnina` works directly. |
| `env` | no | For `.py` scripts only (ignored for binaries). Conda/mamba **env name**: runs via `conda run -n ENV python script.py`. **Path to env directory**: calls that env's `bin/python` directly, no conda needed. |
| `ligand_in_format` | yes | `pdbqt`, `smiles`, or `mol` |
| `ligand_out_format` | yes | `pdbqt`, `sdf`, or `pdb` |
| `input_arg_name` | no | Full CLI flag for the ligand input file (e.g. `--ligand`). If omitted, the ligand string is piped to stdin. |
| `output_arg_name` | no | Full CLI flag for the output file (e.g. `--out`). If omitted, output is read from stdout. |
| `score_mode` | no | `min` (default) or `max` — which end of the score range is "better". |
| `program_args` | no | Dict of extra arguments forwarded to the program as `--key value` pairs. The special key `extra_args` (a plain string) is appended last, split on whitespace — use it for positional arguments or flags that do not follow `--key value` format. |
| `parse_score` | yes | How to extract the docking score (see below). Without it every molecule is silently discarded. |

### `parse_score` sub-keys

| Key | Output format | Description |
|-----|---------------|-------------|
| `score_field` | `sdf` | Name of the SD tag, e.g. `RTMScore`. Supports integer, float, and scientific notation. |
| `score_regex` | `pdbqt` / `pdb` | Regex with one capture group returning the score. Applied to the first MODEL block only. Example: `'REMARK VINA RESULT:\s+(-?[\d.]+)'` |

---

## Example 1 — Vina installed in the active easydock environment

Vina is on the `PATH` of the current conda environment, so no `env` is needed.

```yaml
script_file: vina
ligand_in_format: pdbqt
ligand_out_format: pdbqt
input_arg_name: --ligand
output_arg_name: --out
score_mode: min

program_args:
  receptor: /path/to/protein.pdbqt
  config: /path/to/grid.txt
  exhaustiveness: 8
  num_modes: 5
  cpu: 1
  seed: 0

parse_score:
  score_regex: 'REMARK VINA RESULT:\s+(-?[\d.]+)'
```

---

## Example 2 — CarsiDock installed in a separate conda environment

`env: carsidock_env` tells binary_dock to run the script as
`conda run --no-capture-output -n carsidock_env python /path/to/carsidock_dock.py ...`.
Use a directory path (e.g. `env: ~/miniconda3/envs/carsidock_env`) to bypass conda and call
that environment's Python directly.

```yaml
script_file: /path/to/carsidock_dock.py
env: carsidock_env
ligand_in_format: smiles
ligand_out_format: sdf
input_arg_name: --ligand
output_arg_name: --output
score_mode: max

program_args:
  protein: /path/to/protein.pdb
  reflig: /path/to/reference_ligand.sdf
  num_conformer: 5

parse_score:
  score_field: RTMScore
```

!!! note
    Adjust `program_args` keys and `input_arg_name`/`output_arg_name` to match your specific
    CarsiDock script's CLI.  `env` is ignored if `script_file` does not end in `.py`.

---

## Example 3 — QVina binary file

QVina uses the same argument names as Vina.  Substitute `qvina-w` path for the QVina-W variant.

```yaml
script_file: /path/to/qvina2
ligand_in_format: pdbqt
ligand_out_format: pdbqt
input_arg_name: --ligand
output_arg_name: --out
score_mode: min

program_args:
  receptor: /path/to/protein.pdbqt
  config: /path/to/grid.txt
  exhaustiveness: 8
  num_modes: 5
  cpu: 2
  seed: 0

parse_score:
  score_regex: 'REMARK VINA RESULT:\s+(-?[\d.]+)'
```

---

## Example 4 — Gnina via Apptainer

`script_file` is split on whitespace after path expansion, so the full Apptainer invocation
can be written directly — no wrapper script required.

**Vina scoring** (lower = better):

```yaml
script_file: apptainer run gnina.sif gnina
ligand_in_format: pdbqt
ligand_out_format: pdbqt
input_arg_name: --ligand
output_arg_name: --out
score_mode: min

program_args:
  receptor: /path/to/protein.pdbqt
  config: /path/to/grid.txt
  exhaustiveness: 8
  num_modes: 10
  cpu: 1
  seed: 0
  scoring: vina
  cnn_scoring: none
  addH: 0

parse_score:
  score_regex: 'REMARK VINA RESULT:\s+(-?[\d.]+)'
```

**Default CNN scoring** (higher = better):

```yaml
script_file: apptainer run gnina.sif gnina
ligand_in_format: pdbqt
ligand_out_format: pdbqt
input_arg_name: --ligand
output_arg_name: --out
score_mode: max

program_args:
  receptor: /path/to/protein.pdbqt
  config: /path/to/grid.txt
  exhaustiveness: 8
  num_modes: 10
  cpu: 1
  seed: 0
  scoring: default
  cnn_scoring: rescore
  cnn: dense_ensemble
  addH: 0

parse_score:
  score_regex: 'REMARK CNNscore\s+([\d.]+)'
```
