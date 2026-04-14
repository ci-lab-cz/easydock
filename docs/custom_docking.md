# Adding a Docking Program

This page is intended for **developers** who want to integrate a new docking program into EasyDock by implementing a dedicated `xxx_dock.py` module.

!!! tip "Consider Generic Docking First"
    Before writing a new `xxx_dock.py` module, check whether [Generic Docking](generic_dock.md) meets your needs. Generic docking lets you plug in external binary or Python script through a YAML config, without code changes. Write a dedicated `xxx_dock.py` only when you need custom logic that the generic wrapper cannot express (e.g. batch processing, non-standard I/O, tight integration with a Python library).

## Overview

Each docking program has its own module in the `easydock/` package that exposes a `mol_dock` function:

| Program | Module |
|---|---|
| Vina | `easydock/vina_dock.py` |
| Gnina/Smina | `easydock/gnina_dock.py` |
| QVina | `easydock/qvina_dock.py` |
| Vina-GPU | `easydock/vinagpu_dock.py` |
| Server-based (CarsiDock, SurfDock, Vina-GPU server) | `easydock/server_dock.py` |
| Generic (config-driven wrapper) | `easydock/generic_dock.py` |

## Function Signature

Create a `mol_dock` function with this signature:

```python
def mol_dock(mol, config, ring_sample=False):
    """
    Dock a single molecule.

    Args:
        mol: RDKit Mol object with the _Name property set.
        config: Path to the YAML configuration file.
        ring_sample: If True, dock multiple ring conformers and return the best one.

    Returns:
        tuple: (mol_name, result_dict_or_None)

        result_dict must contain keys matching database field names:
            - docking_score: float — the best-pose score
            - mol_block:     str   — top pose in MDL mol block format
                                     (the first line must be the molecule name)
            - raw_block:     str   — raw output of all poses (optional; PDBQT or SDF)
            - dock_time:     float — wall-clock seconds spent on this molecule (optional)
    """
```

On failure, return `(mol_name, None)`.

## Example Skeleton

```python
import yaml
from rdkit import Chem


def mol_dock(mol, config, ring_sample=False):
    mol_name = mol.GetProp('_Name')

    with open(config) as f:
        config_dict = yaml.safe_load(f)

    try:
        # 1. Convert the RDKit mol to a format your program accepts
        #    (PDBQT, SDF, MOL, SMILES, ...)
        ligand_input = ...

        # 2. Run the docking program (subprocess or native Python)
        ...
        raw_block = ...

        # 3. Parse the best-pose score and structure
        score = ...
        mol_block = ...

        return mol_name, {
            'docking_score': score,
            'mol_block': mol_block,
            'raw_block': raw_block
        }

    except Exception as e:
        print(f'Error docking {mol_name}: {e}')
        return mol_name, None
```

## Configuration File

EasyDock passes the path to a YAML configuration file as the `config` argument. The file layout is up to you; typical contents:

```yaml
script_file: /path/to/my_program
protein: /path/to/protein.pdb
# program-specific parameters
```

## Using the New Module

### Programmatic use

Pass the new `mol_dock` function directly to `docking`:

```python
from easydock.run_dock import docking
from my_module import mol_dock

for mol_id, res in docking(mols,
                           dock_func=mol_dock,
                           dock_config='my_config.yml',
                           ncpu=4):
    ...
```

### Integration with the CLI

To make the new program available via `easydock --program my_program`, add:

1. **A new choice** for `--program` in the argparse definition in `easydock/run_dock.py`.
2. **An import branch** in `run_dock.py` that imports `mol_dock` from your new module when `args.program == 'my_program'`. Search for the existing branches (`if args.program == 'vina': ...`) and add yours alongside.

## Reference Implementations

The existing modules are the best starting point. `vina_dock.py` is the simplest (pure-Python Vina binding); `server_dock.py` and `generic_dock.py` show more elaborate patterns for batched and config-driven execution.