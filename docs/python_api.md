# Python API

EasyDock can be used as a Python library for custom workflows and integration with other tools.

## Basic Docking

### Local Docking

```python
from easydock.run_dock import docking
from easydock.vina_dock import mol_dock
from rdkit import Chem

# Prepare molecules, exactly the same species as provided will be docked 
# (protonate structures preliminary if needed)
smiles = ['CC(=O)O', 'NCC(=O)O', 'NC(C)C(=O)O']
mols = [Chem.MolFromSmiles(smi) for smi in smiles]

# Assign names (required as identifiers)
for mol, smi in zip(mols, smiles):
    mol.SetProp('_Name', smi)

# Run docking
for mol_id, res in docking(mols, dock_func=mol_dock, 
                          dock_config='config.yml', ncpu=4):
    print(f"Molecule: {mol_id}")
    print(f"Results: {res}")
```

Result is a dictionary (keys may vary for different docking programs). Below are the common keys

```python
{'docking_score': -11.856,
 'raw_block': 'string output of all poses in raw docking format (PDBQT or SDF)',
 'mol_block': 'string output of the top pose in MOL format',
 'dock_time': 34}
```

### Using other integrated docking tools

```python
from easydock.gnina_dock import mol_dock
from easydock.qvina_dock import mol_dock
from easydock.vinagpu_dock import mol_dock
from easydock.server_dock import mol_dock
from easydock.generic_dock import mol_dock

# Same code as above, just import other mol_dock function
...
```

### Ring sampling during docking

All `mol_dock` functions have a `ring_sample` argument

```python
# Run docking
for mol_id, res in docking(mols, 
                           dock_func=mol_dock, 
                           dock_config='config.yml',
                           ring_sample=True,
                           ncpu=4):
    print(f"Molecule: {mol_id}")
    print(f"Results: {res}")
```


## Protonation API

There are three groups of implemented protonation functions: file-based, container-based, and Python-based.

### File-based Protonation Functions

These functions use files as input and output. Each program provides two complementary functions: one to run protonation and write the output file, and one to parse that output file.

Currently implemented: Chemaxon (`cxcalc`).

```python
from easydock.protonation import protonate_chemaxon, read_protonate_chemaxon

protonate_chemaxon(input_file_name, output_file_name)
for smi, name in read_protonate_chemaxon(output_file_name):
    print(smi, name)
```

### Container-based Protonation Functions

`protonate_container` launches a container once and streams molecules through it via stdin/stdout. It takes an iterator of `(smiles, name)` tuples and yields protonated `(smiles, name)` tuples. Both Apptainer/Singularity SIF files and Docker images are supported (see [pre-built containers](installation.md#pre-built-containers)).

```python
from easydock.protonation import protonate_container

data = [('NCC(=O)O', 'glycine'), ('C1CNCCN1', 'piperazine')]

# Apptainer/Singularity SIF container
for smi, name in protonate_container(data, container='path/to/unipka.sif'):
    print(smi, name)

# Docker image
for smi, name in protonate_container(data, container='my-protonation-image'):
    print(smi, name)
```

### Python-based Protonation Functions

These functions take an iterable of `(smiles, name)` tuples and return a generator of protonated `(smiles, name)` tuples.

```python
from easydock.protonation import protonate_molgpka, protonate_pkasolver

data = [('NCC(=O)O', 'glycine'), ('C1CNCCN1', 'piperazine')]

for smi, name in protonate_molgpka(data):
    print(smi, name)

for smi, name in protonate_pkasolver(data):
    print(smi, name)
```

!!! warning
If some structures could not be protonated for any reason, they will be omitted in the output. The order of output structures is also not guaranteed.

!!! note
A more detailed description of protonation functions and their requirements is provided in `protonation.py` script 

All those protonate_ family functions can take `pH` argument to specify pH (default 7.4)

```python
for smi, name in protonate_molgpka(data, pH=12):
    print(smi, name)
```
