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

### Implementation Requirements

Implement a protonation function:

```python
def custom_protonate(mol, ph=7.4):
    """
    Protonate a molecule at given pH.
    
    Args:
        mol: RDKit Mol object
        ph: Target pH (default: 7.4)
    
    Returns:
        RDKit Mol object: Protonated molecule
        None: If protonation fails
    """
    pass
```

### Example Implementation

```python
from rdkit import Chem
import subprocess
import tempfile

def my_protonation_tool(mol, ph=7.4):
    """Protonate using external tool."""
    
    # Write input
    with tempfile.NamedTemporaryFile(mode='w', suffix='.smi', delete=False) as f:
        smi = Chem.MolToSmiles(mol)
        f.write(smi)
        input_file = f.name
    
    # Prepare output
    output_file = input_file.replace('.smi', '_out.smi')
    
    # Run protonation tool
    cmd = [
        'my_protonation_tool',
        '--input', input_file,
        '--output', output_file,
        '--ph', str(ph)
    ]
    
    try:
        subprocess.run(cmd, check=True)
        
        # Read result
        with open(output_file) as f:
            protonated_smi = f.readline().strip()
        
        protonated_mol = Chem.MolFromSmiles(protonated_smi)
        
        # Copy properties
        if protonated_mol:
            for prop in mol.GetPropsAsDict():
                protonated_mol.SetProp(prop, mol.GetProp(prop))
        
        return protonated_mol
    
    except Exception as e:
        print(f"Protonation failed: {e}")
        return None
    
    finally:
        # Cleanup
        import os
        os.unlink(input_file)
        if os.path.exists(output_file):
            os.unlink(output_file)
```

### Integration

Add to `protonation.py`:

```python
from .custom_protonation import my_protonation_tool

def protonate_mol(mol, method='molgpka', ph=7.4):
    """Protonate molecule using specified method."""
    
    if method == 'molgpka':
        return molgpka_protonate(mol, ph)
    elif method == 'my_tool':
        return my_protonation_tool(mol, ph)
    # ... other methods
```

## Containerized Tools

### Create Container Function

```python
import subprocess

def run_containerized_tool(container_path, input_file, output_file):
    """
    Run tool in Singularity/Apptainer container.
    
    Args:
        container_path: Path to .sif file
        input_file: Input file path
        output_file: Output file path
    """
    
    cmd = [
        'apptainer', 'exec',
        '-B', f'{os.path.dirname(input_file)}:/data',
        container_path,
        'my_tool',
        '--input', f'/data/{os.path.basename(input_file)}',
        '--output', f'/data/{os.path.basename(output_file)}'
    ]
    
    subprocess.run(cmd, check=True)
```

### Docker on macOS

```python
import platform

def run_container(container_path, *args):
    """Run container with automatic platform detection."""
    
    if platform.system() == 'Darwin':  # macOS
        # Use Docker to run Apptainer container
        cmd = [
            'docker', 'run',
            '--rm',
            '-v', f'{os.getcwd()}:/data',
            'apptainer:latest',
            'exec', container_path
        ] + list(args)
    else:  # Linux
        cmd = ['apptainer', 'exec', container_path] + list(args)
    
    subprocess.run(cmd, check=True)
```

## Database Schema Extension

### Add Custom Fields

Modify database schema to store additional data:

```python
import sqlite3

def extend_database(db_path):
    """Add custom fields to database."""
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Add custom column
    try:
        cursor.execute("""
            ALTER TABLE docking_results 
            ADD COLUMN custom_score REAL
        """)
        conn.commit()
    except sqlite3.OperationalError:
        # Column already exists
        pass
    
    conn.close()
```

### Store Custom Data

```python
def mol_dock_with_custom_data(molecule, config_file):
    """Docking function that returns custom data."""
    
    # ... docking code ...
    
    return mol_name, {
        'docking_score': score,
        'pose_mol': docked_mol,
        'custom_score': custom_score_value  # Custom field
    }
```

## Configuration Parsers

### Custom Config Parser

```python
import yaml

def parse_custom_config(config_file):
    """Parse configuration with custom logic."""
    
    with open(config_file) as f:
        config = yaml.safe_load(f)
    
    # Add computed values
    config['grid_volume'] = (
        config['size_x'] * 
        config['size_y'] * 
        config['size_z']
    )
    
    # Validate
    assert config['exhaustiveness'] > 0, "Invalid exhaustiveness"
    
    return config
```

## Scoring Function Customization

### Custom Scoring

```python
def rescore_poses(mol, pose_mols, protein_file):
    """
    Re-score docked poses with custom function.
    
    Args:
        mol: Original molecule
        pose_mols: List of docked poses
        protein_file: Protein structure file
    
    Returns:
        list: Scores for each pose
    """
    
    scores = []
    for pose in pose_mols:
        # Calculate custom score
        score = calculate_custom_score(pose, protein_file)
        scores.append(score)
    
    return scores
```

## Tips for Customization

### Error Handling

Always include robust error handling:

```python
def safe_mol_dock(molecule, config_file):
    """Docking with error handling."""
    
    mol_name = molecule.GetProp('_Name')
    
    try:
        # Docking code
        result = dock_molecule(molecule, config_file)
        return mol_name, result
    
    except Exception as e:
        print(f"Error docking {mol_name}: {e}", file=sys.stderr)
        return mol_name, None
```

### Logging

Add logging for debugging:

```python
import logging

logger = logging.getLogger(__name__)

def mol_dock_with_logging(molecule, config_file):
    """Docking with detailed logging."""
    
    mol_name = molecule.GetProp('_Name')
    logger.info(f"Starting docking for {mol_name}")
    
    # ... docking code ...
    
    logger.info(f"Finished docking {mol_name}, score: {score}")
    return mol_name, results
```

### Testing

Test custom functions:

```python
def test_custom_dock():
    """Test custom docking function."""
    
    from rdkit import Chem
    
    mol = Chem.MolFromSmiles('CCO')
    mol.SetProp('_Name', 'test')
    
    mol_name, results = my_custom_dock(mol, 'test_config.yml')
    
    assert mol_name == 'test'
    assert 'docking_score' in results
    assert results['pose_mol'] is not None
```
