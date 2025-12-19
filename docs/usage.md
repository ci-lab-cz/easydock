# Usage Guide

There are two major stages which can run separately or simultaneously: ligand preparation and their docking.

## Database Initialization

EasyDock uses an SQLite database to store all inputs, parameters, and results. Once initialized, the database preserves all settings for reproducibility.

!!! info "Database Immutability"
    Once a database is created, most parameters cannot be changed via command line. The database is never overwritten. Delete and reinitialize if settings are incorrect. Direct changes of settings in the database should be done carefully to not create inconsistent settings.

### Basic Initialization

Create a database with molecule validation:

```bash
easydock -i input.smi -o output.db -c 4
```

**Parameters:**

- `-i`: Input SMILES file
- `-o`: Output database file
- `-c`: Number of CPU cores

This performs:

- Conversion to 3D structures (or uses existing coordoinates if 3D SDF input)
- Salt removal
- Validation

### Stereoisomer Enumeration

Enumerate up to 4 stereoisomers for undefined chiral centers and double bonds:

```bash
easydock -i input.smi -o output.db -c 4 -s 4
```

Default is 1 stereoisomer (reproducible enumeration). Maximum runtime per molecule is limited to 300 seconds. This avoids excessively long generation for some structures. In those cases it is recommended to enumerate stereoisomers before supplying these molecules to EasyDock.

### Protonation Options

**Using MolGpKa:**
```bash
easydock -i input.smi -o output.db -c 4 --protonation molgpka
```

**Using Uni-pKa container:**
```bash
easydock -i input.smi -o output.db -c 4 --protonation /path/to/unipka.sif
```

**Using Chemaxon (requires license):**
```bash
easydock -i input.smi -o output.db -c 4 --protonation chemaxon
```

**No protonation (use input states):**
```bash
easydock -i input.smi -o output.db -c 4
```

## Molecular Docking

### Vina Docking

**Complete pipeline (initialization + docking):**
```bash
easydock -i input.smi -o output.db --program vina --config config.yml -c 4 --sdf
```

**Using pre-initialized database:**
```bash
easydock -o output.db --program vina --config config.yml -c 4 --sdf
```

**Vina Configuration (config.yml):**
```yaml
protein: /path/to/protein.pdbqt
protein_setup: /path/to/grid.txt
exhaustiveness: 8
seed: 0
n_poses: 5
ncpu: 5
```

**Grid Box Definition (grid.txt):**
```
center_x = 10.0
center_y = 15.0
center_z = 20.0
size_x = 25.0
size_y = 25.0
size_z = 25.0
```

!!! note "CPU Parameters"
    - Command line `-c 4`: Docks 4 molecules in parallel
    - Config `ncpu: 5`: Uses 5 CPUs per molecule
    - Total CPUs: 4 × 5 = 20
    
    Ensure the product matches or slightly exceeds available CPUs.

### Gnina Docking

```bash
easydock -i input.smi -o output.db --program gnina --config config.yml -c 4 --sdf
```

**Gnina Configuration:**
```yaml
script_file: /path/to/gnina
protein: /path/to/protein.pdbqt
protein_setup: /path/to/grid.txt
exhaustiveness: 8
scoring: default
cnn_scoring: rescore
cnn: dense_ensemble
n_poses: 10
addH: False
ncpu: 1
seed: 0
```

### Smina Docking

Use Gnina program with Smina-specific configuration:

```yaml
script_file: /path/to/gnina
protein: /path/to/protein.pdbqt
protein_setup: /path/to/grid.txt
exhaustiveness: 8
scoring: vinardo
cnn_scoring: None
cnn: dense_ensemble
n_poses: 10
addH: False
ncpu: 1
seed: 0
```

### Vina-GPU Family

For Vina-GPU, QVina2-GPU, or QVinaW-GPU:

```bash
easydock -i input.smi -o output.db --program vina-gpu --config config.yml --sdf
```

**Configuration:**
```yaml
script_file: /path/to/AutoDock-Vina-GPU-2-1 --opencl_binary_path /path/to/opencl/
protein: /path/to/protein.pdbqt
protein_setup: /path/to/grid.txt
n_poses: 3
thread: 8000  # Optional, default is 8000
```

!!! tip "GPU Program Variants"
    Change `script_file` to use different variants:
    
    - `AutoDock-Vina-GPU-2-1` for Vina-GPU
    - `QuickVina2-GPU-2-1` for QVina2-GPU
    - `QuickVina-W-GPU-2-1` for QVinaW-GPU

### QVina2 / QVina-W

```bash
easydock -i input.smi -o output.db --program qvina --config config.yml --sdf
```

**Configuration:**
```yaml
script_file: /path/to/qvina2  # or /path/to/qvinaw
protein: /path/to/protein.pdbqt
protein_setup: /path/to/grid.txt
exhaustiveness: 8
n_poses: 3
ncpu: 2
seed: 0
```

## Resuming Interrupted Calculations

Simply run the same command or provide just the database:

```bash
easydock -o output.db
```

EasyDock will continue from where it stopped, using settings stored in the database.

## Output Options

### Generate SDF File

Add `--sdf` flag to create an SDF file with top poses:

```bash
easydock -o output.db --program vina --config config.yml -c 4 --sdf
```

This creates `output.sdf` with the best scoring pose for each molecule.
