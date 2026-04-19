# Generic Docking

`--program generic` is a generic wrapper that can run **external docking binary or Python script** through a single YAML config file. No code changes are needed — the config drives everything: which executable to call, how to pass the ligand and receive the output, and how to extract the score.

```bash
easydock -i input.smi -o output.db --program generic --config config.yml -c 4 --sdf
```

The main conventions:
1. A program or script takes a single molecule to dock at a single call
2. A program or script returns a file with docked poses in PDBQT/PDB/SDF format ordered by the score (the top scored pose in the first one)

## Config reference

| Key | Required | Description                                                                                                                                                                                                                                      |
|-----|----------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `script_file` | yes | Executable or command string. `~` and environment variables are <br>expanded first, then the result is split on whitespace, <br>so `apptainer run gnina.sif gnina` works directly.                                                                   |
| `env` | no | For `.py` scripts only (ignored for binaries). Conda/mamba **env name**: <br>runs via `conda run -n ENV python script.py`. <br>**Path to env directory**: calls that env's `bin/python` directly, <br>no conda needed.                                       |
| `ligand_in_format` | yes | `pdbqt`, `smiles`, or `mol`                                                                                                                                                                                                                      |
| `ligand_out_format` | yes | `pdbqt`, `sdf`, or `pdb`                                                                                                                                                                                                                         |
| `input_arg_name` | no | Full CLI flag for the ligand input file (e.g. `--ligand`). <br>If omitted, the ligand string is piped to stdin.                                                                                                                                      |
| `output_arg_name` | no | Full CLI flag for the output file (e.g. `--out`). <br>If omitted, output is read from stdout.                                                                                                                                                        |
| `score_mode` | no | `min` (default) or `max` — which end of the score range is "better".                                                                                                                                                                             |
| `program_args` | no | Dict of extra arguments forwarded to the program as `--key value` pairs. <br>The special key `extra_args` (a plain string) is appended last, split on whitespace — <br>use it for positional arguments or flags that do not follow <br>`--key value` format. |
| `parse_score` | yes | How to extract the docking score (see below). Without it every <br>molecule is silently discarded.                                                                                                                                                   |

### `parse_score` sub-keys

| Key | Output format | Description |
|-----|---------------|-------------|
| `score_field` | `sdf` | Name of the SD tag, e.g. `RTMScore`. Supports integer, float, and <br>scientific notation. |
| `score_regex` | `pdbqt` / `pdb` | Regex with one capture group returning the score. Applied to the first <br>MODEL block only. Example: `'REMARK VINA RESULT:\s+(-?[\d.]+)'` |

---

## Example 1 — QVina binary file

QVina uses input and output file format PDBQT. Those files are supplied with `--input` and `--out` arguments. The lower the docking score the better (`score_mode: min`). 
`program_args` further arguments which will be passed to the program "as is". Regex `REMARK VINA RESULT:\s+(-?[\d.]+)` will extract the relevant docking score value from the first item of the output file.   

```yaml
script_file: ~/qvina/bin/qvina2.1
ligand_in_format: pdbqt
ligand_out_format: pdbqt
input_arg_name: --ligand
output_arg_name: --out
score_mode: min

program_args:
  receptor: ~/2btr/A/2btr_A.pdbqt
  config: ~/2btr/grid.txt
  exhaustiveness: 8
  num_modes: 5
  cpu: 2
  seed: 0

parse_score:
  score_regex: REMARK VINA RESULT:\s+(-?[\d.]+)
```

---

## Example 2 (generic) — a program run via Apptainer

A user can provide a full command line call. `script_file` is split on whitespace after path expansion, so the full Apptainer invocation
can be written directly — no wrapper script required.

```yaml
script_file: apptainer run container.sif program_name
ligand_in_format: smiles
ligand_out_format: pdbqt
input_arg_name: --ligand
output_arg_name: --out
score_mode: min

program_args:
  receptor: /path/to/protein.pdbqt
  config: /path/to/grid.txt
  exhaustiveness: 8
  seed: 10

parse_score:
  score_regex: 'REMARK VINA RESULT:\s+(-?[\d.]+)'
```

---

## Example 3 (generic) — a program installed in a separate conda environment

`env: program_env` tells generic_dock to run the script as 
`conda run --no-capture-output -n program_env python /path/to/docking_program.py ...`. 
Use a directory path (e.g. `env: ~/miniconda3/envs/program_env`) 
to bypass conda and call that environment's Python directly.
Since output is in SDF format, the docking score will be extracted from the field `RTMScore`.   

```yaml
script_file: /path/to/docking_program.py
env: program_env
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
    `env` is ignored if `script_file` does not end in `.py`. If `script_file` ends with `.py` and no `env` was provided, the script will be executed in a current Python environment. 
