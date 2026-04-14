# Adding a Protonation Tool

This page is intended for **developers** who want to integrate a new protonation tool into EasyDock. Three conventions are supported — pick the one that matches how the underlying tool is distributed:

| Convention | When to use | Examples |
|---|---|---|
| [**File-based**](#file-based-protonation) | Tool is an external binary that reads/writes files | Chemaxon (`cxcalc`) |
| [**Native Python**](#native-python-protonation) | Tool is a Python library; pure-Python workflow | MolGpKa, pkasolver |
| [**Container-based**](#container-based-protonation) | Tool has complex dependencies; ship it in Apptainer/Docker | Uni-pKa |

All integrations are centralized in `easydock/protonation.py`. Registration with the CLI happens in the `add_protonation` function in `easydock/database.py`.

---

## File-based Protonation

Use this convention when the tool is a command-line program that takes an input file and writes an output file.

### Two functions required

1. **`protonate_xxx(input_fname, output_fname, ...)`** — runs the external program.
    - Input: tab-separated SMILES file (columns: `smiles`, `mol_name`).
    - Output: any format you can read back in step 2.

2. **`read_protonate_xxx(fname)`** — generator yielding `(smiles, mol_name)` tuples parsed from the output file.

### Example

Reference implementation: `protonate_chemaxon` / `read_protonate_chemaxon` in `easydock/protonation.py`.

```python
import subprocess
from rdkit import Chem


def protonate_xxx(input_fname, output_fname, pH: float = 7.4):
    subprocess.run(['my_tool', '--pH', str(pH),
                    '-i', input_fname, '-o', output_fname], check=True)


def read_protonate_xxx(fname):
    # parse the output file produced by my_tool
    for line in open(fname):
        smi, name = line.strip().split('\t')
        yield smi, name
```

### Registration in `add_protonation`

In `easydock/database.py`, add a branch in `add_protonation` matching the existing `'chemaxon'` pattern:

```python
elif program == 'my_tool':
    protonate_func = partial(protonate_xxx, pH=pH)
    read_func = read_protonate_xxx
    # (chunked file loop, same as for chemaxon)
```

---

## Native Python Protonation

Use this convention when the tool is a pure-Python library (RDKit-based, PyTorch model, etc.). Avoid file I/O — work directly on `(smi, name)` tuples.

### Function signature

A single generator:

```python
def protonate_xxx(items: Iterator[Tuple[str, str]],
                  ncpu: int = 1,
                  pH: float = 7.4) -> Iterator[Tuple[str, str]]:
    """Take (smi, name) tuples, yield (protonated_smi, name) tuples."""
```

- Use `multiprocessing.Pool` internally when the work parallelizes cleanly (see `protonate_pkasolver`).
- If parallel execution is slower than single-process (as with MolGpKa), just iterate.
- On error, log a warning and yield `(None, name)` — the downstream database layer will skip that molecule.

### Example

Reference implementations: `protonate_molgpka`, `protonate_pkasolver` in `easydock/protonation.py`.

```python
def protonate_xxx(items, ncpu: int = 1, pH: float = 7.4):
    from my_tool import predict_major_microspecies

    for smi, name in items:
        try:
            new_smi = predict_major_microspecies(smi, pH=pH)
            yield new_smi, name
        except Exception:
            logging.warning(f'{name} caused an error during protonation, skipping')
            yield None, name
```

### Registration in `add_protonation`

Add a branch matching the existing `'molgpka'` / `'pkasolver'` pattern in `easydock/database.py`:

```python
elif program == 'my_tool':
    protonate_func = partial(protonate_xxx, ncpu=ncpu, pH=pH)
    # (streaming loop, same as for molgpka)
```

Also add the program name to `protonation_programs` in `easydock/args_validation.py` so the CLI accepts it.

---

## Container-based Protonation

Use this convention when the tool has heavy or conflicting dependencies (specific CUDA, conda environment, etc.) and is easier to ship as a container. Both **Apptainer/Singularity SIF** files and **Docker images** are supported through the same interface.

EasyDock already provides a generic driver — `protonate_container` in `easydock/protonation.py` — so you usually do **not** need to write Python glue code. You only need to build a container that conforms to the interface below.

### Container interface

The container must expose a `protonate` subcommand that:

1. Reads lines from **STDIN**, each line formatted as `smiles\tname\n`.
2. For each input molecule, writes the protonated result to **STDOUT** as `protonated_smiles\tname\n`.
3. Accepts a `--pH <float>` argument (required; default convention 7.4).
4. Keeps running until STDIN is closed — the container is launched **once** and streams all molecules through a single process.

EasyDock invokes the container as:

```bash
# Apptainer / Singularity
apptainer run [--nv] /path/to/container.sif protonate --pH 7.4

# Docker
docker run -i [--gpus all] my-image protonate --pH 7.4
```

`--nv` / `--gpus all` is added automatically when `nvidia-smi` is available. `-i` is always added for Docker so STDIN stays open.

### Apptainer `%runscript` template

```singularity
%runscript
    case "$1" in
        protonate)
            shift
            exec python /opt/mytool/protonate.py "$@"
            ;;
        help)
            echo "Usage: apptainer run mytool.sif protonate --pH 7.4"
            ;;
        *)
            exec "$@"
            ;;
    esac
```

### Docker `ENTRYPOINT` template

```dockerfile
ENTRYPOINT ["/bin/sh", "-c", "\
  case \"$1\" in \
    protonate) shift; exec python /opt/mytool/protonate.py \"$@\" ;; \
    help)      echo 'Use: docker run <image> protonate --pH 7.4' ;; \
    *)         exec \"$@\" ;; \
  esac", "--"]
```

### Inner script skeleton

The `protonate.py` inside the container just reads STDIN and writes STDOUT — no files, no argparse for `-i`/`-o` required:

```python
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--pH', type=float, default=7.4)
args = parser.parse_args()

for line in sys.stdin:
    smi, name = line.rstrip('\n').split('\t')
    new_smi = my_protonation(smi, pH=args.pH)
    sys.stdout.write(f'{new_smi}\t{name}\n')
    sys.stdout.flush()   # important: flush after every line
```

!!! warning "Flush every line"
    EasyDock reads one line at a time from the container's STDOUT. Without `sys.stdout.flush()` after each molecule, Python's default block-buffered stdout will stall the pipeline.

### Reference implementation

See `containers/unipka/` in the repository for a complete example (`unipka.def`, `Dockerfile`, `unipka.py`).

### Using the container

No code changes are needed in EasyDock. Users pass the SIF path or Docker image name directly via `--protonation`:

```bash
# SIF container
easydock -i input.smi -o output.db -c 4 --protonation /path/to/mytool.sif

# Docker image
easydock -i input.smi -o output.db -c 4 --protonation my-protonation-image
```

The dispatcher in `add_protonation` (`easydock/database.py`) detects SIF files by extension and Docker images by elimination (any `--protonation` value that is neither a known built-in name nor an existing non-SIF file is treated as a Docker image name) and routes the call through `protonate_container`.