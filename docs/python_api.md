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
from easydock.qnina_dock import mol_dock
from easydock.vinagpu_dock import mol_dock

# Same code as above, just import other mol_dock function
...
```

## Ring sampling

All `mol_dock` functions have a `ring_sample` argument

```python
import partial

# Run docking
for mol_id, res in docking(mols, 
                           dock_func=partial(mol_dock, ring_sample=True), 
                           dock_config='config.yml', ncpu=4):
    print(f"Molecule: {mol_id}")
    print(f"Results: {res}")
```










## Working with Results

### Access Result Dictionary

```python
for mol_id, res in docking(mols, dock_func=mol_dock,
                          dock_config='config.yml', ncpu=4):
    if res is not None:
        print(f"ID: {mol_id}")
        print(f"Score: {res.get('docking_score')}")
        print(f"Pose: {res.get('pose_mol')}")  # RDKit Mol object
```

### Save Results

```python
from rdkit import Chem

results = []
for mol_id, res in docking(mols, dock_func=mol_dock,
                          dock_config='config.yml', ncpu=4):
    if res is not None:
        results.append((mol_id, res))

# Save to SDF
writer = Chem.SDWriter('output.sdf')
for mol_id, res in results:
    pose = res['pose_mol']
    pose.SetProp('docking_score', str(res['docking_score']))
    pose.SetProp('mol_id', mol_id)
    writer.write(pose)
writer.close()
```

## Database Integration

### Initialize Database

```python
from easydock.database import initialize_db
from rdkit import Chem

# Create molecules
mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
for i, mol in enumerate(mols):
    mol.SetProp('_Name', f'mol_{i}')

# Initialize database
initialize_db(
    mols=mols,
    output_db='output.db',
    ncpu=4,
    protonation='molgpka',
    max_stereoisomers=2
)
```

### Query Database

```python
import sqlite3

conn = sqlite3.connect('output.db')
cursor = conn.cursor()

# Get top scoring molecules
cursor.execute("""
    SELECT mol_id, docking_score 
    FROM docking_results 
    WHERE pose_num = 1 
    ORDER BY docking_score ASC 
    LIMIT 10
""")

for mol_id, score in cursor.fetchall():
    print(f"{mol_id}: {score}")

conn.close()
```

## Custom Workflows

### Pre-filter Molecules

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

# Filter by molecular weight
filtered_mols = []
for mol in mols:
    if 200 < Descriptors.MolWt(mol) < 500:
        filtered_mols.append(mol)

# Dock filtered molecules
results = docking(filtered_mols, dock_func=mol_dock,
                 dock_config='config.yml', ncpu=4)
```

### Post-process Results

```python
from rdkit.Chem import AllChem

good_poses = []
for mol_id, res in docking(mols, dock_func=mol_dock,
                          dock_config='config.yml', ncpu=4):
    if res and res['docking_score'] < -8.0:
        pose = res['pose_mol']
        # Calculate RMSD or other properties
        good_poses.append((mol_id, pose, res['docking_score']))

# Sort by score
good_poses.sort(key=lambda x: x[2])
```

## Custom Docking Function

Implement your own docking function:

```python
def custom_mol_dock(molecule, config_file):
    """
    Custom docking function.
    
    Args:
        molecule: RDKit Mol object with _Name property
        config_file: Path to YAML configuration
    
    Returns:
        tuple: (molecule_name, results_dict)
    """
    import yaml
    from rdkit import Chem
    
    # Get molecule name
    mol_name = molecule.GetProp('_Name')
    
    # Load config
    with open(config_file) as f:
        config = yaml.safe_load(f)
    
    # Your docking code here
    # ...
    
    # Return results
    results = {
        'docking_score': score,
        'pose_mol': docked_mol,
        'custom_field': custom_value
    }
    
    return mol_name, results

# Use custom function
for mol_id, res in docking(mols, dock_func=custom_mol_dock,
                          dock_config='config.yml', ncpu=4):
    print(mol_id, res)
```

## Protonation API

### Use Protonation Functions

```python
from easydock.protonation import protonate_mol

mol = Chem.MolFromSmiles('CC(=O)O')
protonated = protonate_mol(mol, method='molgpka')
```

### Batch Protonation

```python
from easydock.protonation import batch_protonate

mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
protonated_mols = batch_protonate(mols, method='molgpka', ncpu=4)
```

## Error Handling

```python
from easydock.run_dock import docking
from easydock.vina_dock import mol_dock

successful = []
failed = []

for mol_id, res in docking(mols, dock_func=mol_dock,
                          dock_config='config.yml', ncpu=4):
    if res is not None:
        successful.append((mol_id, res))
    else:
        failed.append(mol_id)
        print(f"Docking failed for {mol_id}")

print(f"Successful: {len(successful)}")
print(f"Failed: {len(failed)}")
```

## Advanced: Distributed Docking

```python
from dask.distributed import Client
from easydock.run_dock import docking_distributed
from easydock.vina_dock import mol_dock

# Connect to Dask cluster
client = Client('scheduler-address:8786')

# Run distributed docking
results = docking_distributed(
    mols=mols,
    dock_func=mol_dock,
    dock_config='config.yml',
    client=client
)

for mol_id, res in results:
    print(mol_id, res)

client.close()
```

---
