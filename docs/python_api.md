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
 'pdb_block': 'string output of all poses in PDB format',
 'mol_block': 'string output of the top pose in MOL format',
 'dock_time': 34}
```

### Using other integrated docking tools

```python
from easydock.gnina_dock import mol_dock
from easydock.qvina_dock import mol_dock
from easydock.vinagpu_dock import mol_dock
from easydock.server_dock import mol_dock

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

There are two groups of implemented protonation functions: file-based and Python-based. The former use files as input and output, the latter use Python objects.

### File-based Protonation Functions

Currently, these are functions to protonate using Chemaxon or sif-containers, e.g. uni-pka.sif ([list](installation.md#pre-build-containers)).
Input: tab-separated SMILES file with molecule names.
Output: tab-separated SMILES file for major microspecies with molecule names.
Each program has two complementary functions: to protonate structure in input file and return output file, and a function to parse the output file to return. 

```python
from easydock.protonation import protonate_chemaxon, protonate_apptainer, read_protonate_chemaxon, read_smiles

protonate_chemaxon(input_file_name, output_file_name)
for smi, name in read_protonate_chemaxon(output_file_name):
    print(smi, name)

protonate_apptainer(input_file_name, output_file_name, 'path/to/container.sif')
for smi, name in read_smiles(output_file_name):
    print(smi, name)
```

### Python-based Protonation Functions

These functions take an iterable of tuples containing SMILES and a molecule name as input and return a generator over SMILES of a predicted major microspecies and a name

```python
from easydock.protonation import protonate_molgpka, protonate_pkasolver

data = [('NCC(=O)O', 'glycine'), ('C1CNCCN1', 'piperazine')]

for smi. name in protonate_molgpka(data):
    print(smi, name)

for smi. name in protonate_pkasolver(data):
    print(smi, name)
```

!!! warning
If some structures could not be protonated for any reason, they will be omitted in the output. The order of output structures is also not guaranteed.

!!! note
A more detailed description of protonation functions and their requirements is provided in `protonation.py` script 

All those protonate_ family functions can take `pH` argument to specify pH (default 7.4)

```python
for smi. name in protonate_molgpka(data, pH=12):
    print(smi, name)
```
