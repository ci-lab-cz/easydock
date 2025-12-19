# Python API

EasyDock can be used as a Python library for custom workflows and integration with other tools.

## Basic Docking

### Local Docking

```python
from easydock.run_dock import docking
from easydock.vina_dock import mol_dock
from rdkit import Chem

# Prepare molecules
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

### Use Protonation Functions

### Batch Protonation
