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

## Reusable Persistent Server Backend

EasyDock includes a reusable backend (`--program server`) for tools with expensive model initialization.
It starts one long-lived server process per worker and reuses it for multiple molecules.

### When to use

- Containerized tools launched with Apptainer/Singularity
- ML-based docking where loading model weights dominates per-molecule runtime
- Backends that can expose request/response API over stdin/stdout

### Required protocol (JSON Lines)

Each request and response is one JSON object per line:

```json
{"id": 1, "command": "init", "payload": {"protein": "/path/protein.pdbqt"}}
{"id": 1, "status": "ok"}
{"id": 2, "command": "dock_batch", "payload": {"molecule_id": "mol_1", "ligands_pdbqt": ["..."]}}
{"id": 2, "status": "ok", "results": [{"docking_score": -8.3, "pdb_block": "MODEL ... ENDMDL"}]}
```

### Config keys

- `script_file`: command used to start server (supports complex shell commands, including `apptainer exec ...`)
- `init_command`: init command name (default `init`)
- `dock_command`: docking command name (default `dock_batch`)
- `result_items_key`: key containing list of conformer-level outputs (default `results`)
- `score_key`: key with numeric score in each result item (default `docking_score`)
- `pose_key`: key with PDBQT pose text in each result item (default `pdb_block`)
- `mol_block_key`: key with ready-to-store MOL block in each result item (default `mol_block`)
- `ligand_payload_type`: payload format sent to server (`pdbqt`, `smiles`, `mol_block`; default `pdbqt`)
- `ligand_payload_key`: payload key for ligand list (default depends on `ligand_payload_type`)
- `molecule_id_payload_key`: key for molecule id in dock request payload (default `molecule_id`)
- `ring_sample_payload_key`: key for ring sampling flag in dock request payload (default `ring_sample`)
- `score_mode`: `min` (default) or `max`
- `startup_timeout`: server startup timeout in seconds
- `request_timeout`: per-request timeout in seconds
- `init_payload`: optional dict sent during init; if omitted, all non-control config values are sent automatically

## Custom Protonation Tools

### Containerized solutions

EasyDock supports singularity/apptainer containers (SIF-files) equipped with protonation tools and all dependencies. The container should provide a command "protonate" which takes two necessary arguments -i/--input and -o/--output taking input and output files in SMILES format: two columns tab-separated with SMILES and molecule id.
Therefore, the protonation should be invoked by a command:
```bash
apptainer run container.sif protonate -i input.smi -o output.smi 
```
