# Usage Guide

There are two major stages which can run separately or simultaneously: ligand preparation and their docking.

## Database Initialization / Ligand Preparation

EasyDock uses an SQLite database to store all inputs, parameters, and results. Once initialized, the database preserves all settings for reproducibility.

!!! info "Database Immutability"
    Once a database is created, most parameters cannot be changed via command line. The database is never overwritten. Delete and reinitialize if settings are incorrect. Direct changes of settings in the database (`setup` table) should be done carefully to not create inconsistent settings.

### Basic Initialization

Create a database with molecule validation:

```bash
easydock -i input.smi -o output.db -c 4
```

Parameters:
- `-i`: Input SMILES file  
- `-o`: Output database file  
- `-c`: Number of CPU cores  

This performs:
- Salt removal  
- Generation of one stereoisomer if there are some undefined chiral centers or double bonds (reproducible generation)  
- Conversion to 3D structures (or uses existing coordinates if 3D SDF input)  

### Stereoisomer Enumeration

Enumerate up to 4 stereoisomers for undefined chiral centers and double bonds:

```bash
easydock -i input.smi -o output.db -c 4 -s 4
```

Maximum runtime of stereoisomer enumeration is limited to 300 seconds per molecule. This avoids excessively long generation for some structures. In those cases it is recommended to enumerate stereoisomers before supplying these molecules to EasyDock.

### Protonation Options

Using MolGpKa:
```bash
easydock -i input.smi -o output.db -c 4 --protonation molgpka
```

Using pre-built Uni-pKa container:
```bash
easydock -i input.smi -o output.db -c 4 --protonation /path/to/unipka.sif
```

Using Chemaxon (requires license):
```bash
easydock -i input.smi -o output.db -c 4 --protonation chemaxon
```

Custom pH
```bash
easydock -i input.smi -o output.db -c 4 --protonation molgpka --pH 12
```

No protonation (use input states):
```bash
easydock -i input.smi -o output.db -c 4
```

## Molecular Docking

### Vina Docking

Complete pipeline (initialization + docking):
```bash
easydock -i input.smi -o output.db --program vina --config config.yml -c 4 --sdf
```

Using pre-initialized database:
```bash
easydock -o output.db --program vina --config config.yml -c 4 --sdf
```

Vina Configuration (config.yml):
```yaml
protein: /path/to/protein.pdbqt
protein_setup: /path/to/grid.txt
exhaustiveness: 8
seed: 0
n_poses: 5
ncpu: 5
```

Grid Box Definition (grid.txt):
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

Gnina Configuration:
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

Configuration:
```yaml
script_file: /path/to/AutoDock-Vina-GPU-2-1 --opencl_binary_path /path/to/opencl/
protein: /path/to/protein.pdbqt
protein_setup: /path/to/grid.txt
n_poses: 3
thread: 8000  # Optional, default is 8000
```

!!! tip "Program Variants"
    Change `script_file` to use different variants:
    
    - `AutoDock-Vina-GPU-2-1` for Vina-GPU
    - `QuickVina2-GPU-2-1` for QVina2-GPU
    - `QuickVina-W-GPU-2-1` for QVinaW-GPU

### QVina2 / QVina-W

```bash
easydock -i input.smi -o output.db --program qvina --config config.yml --sdf
```

Configuration:
```yaml
script_file: /path/to/qvina2  # or /path/to/qvinaw
protein: /path/to/protein.pdbqt
protein_setup: /path/to/grid.txt
exhaustiveness: 8
n_poses: 3
ncpu: 2
seed: 0
```

!!! tip "Program Variants"
    Change `script_file` to use different variants:
    
    - `qvina2`
    - `qvina-w`
    - `vina` - Vina can be used through a binary executable

### Server-Based Docking

`--program server` enables docking through a long-running server process. The server is started once per worker and handles multiple molecules without restarting. This mode supports containerized docking programs such as CarsiDock and GPU-accelerated Vina variants.

```bash
easydock -i input.smi -o output.db --program server --config config.yml -c 4 --sdf
```

!!! tip "GPU Auto-Detection"
    If `nvidia-smi` is accessible, EasyDock automatically adds `--nv` (Apptainer/Singularity) or `--gpus all` (Docker) to the container launch command. No manual configuration is required.

The config file has three kinds of keys:

- **Control keys** (`script_file`, `score_mode`, `startup_timeout`, etc.) — read by EasyDock itself and placed at the top level.
- **`init_server:`** — parameters forwarded verbatim to the server's `init` command (receptor files, pocket definition, program-specific settings).
- **`info_server:`** — optional overrides for the server's default INFO values such as `batch_size`. The server exposes these defaults via its `info` command; this section lets you change them without modifying the server.

The ligand input/output formats (`ligand_in_format`, `ligand_out_format`) are auto-detected from the server's `info` command and do not need to be set manually.

#### CarsiDock

CarsiDock is a deep-learning docking program. It accepts SMILES input and uses RTMScore for pose ranking (higher score = better). A pre-built Apptainer container is available (see [Installation](installation.md#pre-build-containers)).

```yaml
script_file: /path/to/carsidock.sif
score_mode: max   # RTMScore: higher score is better

init_server:
  protein: /path/to/protein.pdb
  reflig: /path/to/reference_ligand.sdf
  num_conformer: 5
```

- `protein`: protein PDB file (defines the binding site)
- `reflig`: reference ligand that defines the binding pocket
- `num_conformer`: number of conformers to generate per ligand (default: 5)

To limit the number of molecules docked per server request (e.g. for memory reasons), override `batch_size` via `info_server`:

```yaml
script_file: /path/to/carsidock.sif
score_mode: max

init_server:
  protein: /path/to/protein.pdb
  reflig: /path/to/reference_ligand.sdf
  num_conformer: 5

info_server:
  batch_size: 5
```

#### Vina-GPU Server

The Vina-GPU server bundles Vina-GPU, QVina2-GPU, and QVinaW-GPU in a single container. It accepts PDBQT input and output and uses Vina scoring (lower = better).

```yaml
script_file: /path/to/vinagpu.sif

init_server:
  protein: /path/to/protein.pdbqt
  protein_setup: /path/to/grid.txt
  program: vina-gpu
  n_poses: 9
  thread: 8000
  seed: 0
```

!!! tip "Program Variants"
    Change the `program` field under `init_server` to switch between GPU variants:

    - `vina-gpu` — AutoDock Vina-GPU
    - `qvina-gpu` — QuickVina2-GPU
    - `qvinaw-gpu` — QuickVina-W-GPU

## Resuming Interrupted Calculations

Simply run the same command or provide just the database:

```bash
easydock -o output.db
```

EasyDock will continue from where it stopped, using settings stored in the database.

## Output Options

### Generate SDF File

Add `--sdf` flag to create an SDF file with top poses after docking will be finished:

```bash
easydock -o output.db --program vina --config config.yml -c 4 --sdf
```

This creates `output.sdf` with the best scoring pose for each molecule.

!!! note "Feature"
    The argument `--sdf` automatically extracts only one stereoisomer with the best docking scores among generated ones by EasyDock. If different stereoisomers were supplied as input, they will be treated as individual species by `--sdf` option. 
