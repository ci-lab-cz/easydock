# Customization

EasyDock can be extended to support custom docking programs and protonation tools.

## Custom Docking Programs

### Implementation Requirements

Create a `mol_dock` function with this signature:

```python
def mol_dock(mol, config, ring_sample=False):
    """
    Dock a single molecule.
    
    Args:
        mol: RDKit Mol object with _Name property set
        config: Path to YAML configuration file
    
    Returns:
        tuple: (molecule_name, results_dict)
        
        results_dict must contain keys matching database field names:
        - docking_score: float
        - mol_block: string representation of the top pose in MOL format
        - other fields as needed
    """
    pass
```

### Example Implementation

```python
import yaml
from rdkit import Chem

def my_docking_program_dock(mol, config):

    # Get molecule identifier
    mol_name = mol.GetProp('_Name')
    
    # Load configuration
    with open(config) as f:
        config_dict = yaml.safe_load(f)
    
    # Prepare ligand file
    # convert a ligand to a suitable format
    ...
    
    try:
        # Run docking through subprocess or natively from python
        ...
        
        # Parse output
        ...
        
        # Return results
        return mol_name, {
            'docking_score': score,
            'mol_block': mol_block
        }
    
    except Exception as e:
        print(f"Error docking {mol_name}: {e}")
        return mol_name, None
```

### Configuration File

Create corresponding YAML configuration:

```yaml
script_file: /path/to/my_docking_program
protein: /path/to/protein.pdb
# Add program-specific parameters
```

### Usage

```python
from easydock.run_dock import docking
from my_module import my_docking_program_dock

results = docking(
    mols=molecules,
    dock_func=my_docking_program_dock,
    dock_config='my_config.yml',
    ncpu=4
)
```

## Custom Protonation Tools

### Containerized solutions

EasyDock supports singularity/apptainer containers (SIF-files) equipped with protonation tools and all dependencies. The container should provide a command "protonate" which takes two necessary arguments -i/--input and -o/--output taking input and output files in SMILES format: two columns tab-separated with SMILES and molecule id.
Therefore, the protonation should be invoked by a command:
```bash
apptainer run container.sif protonate -i input.smi -o output.smi 
```
