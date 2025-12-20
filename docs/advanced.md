# Advanced Features

## Ring Conformer Sampling

Improve docking of compounds containing saturated rings by sampling multiple starting conformers.

```bash
easydock -i input.smi -o output.db --program vina --config config.yml -c 4 --ring_sample
```

How it works:

1. Generate multiple representative conformers for saturated rings
2. Dock each initial conformer independently
3. Store only the best scoring conformer

!!! warning "Performance Impact"
    - Increases computational time proportionally to the number of rings
    - May result in much longer simulation times for compounds with multiple saturated rings

## Containerized Programs

Run docking programs through Docker/Singularity/Apptainer containers.

Example with Gnina in Apptainer:

```yaml
script_file: apptainer exec -B /path/to/data:/path/to/data container.sif gnina
protein: /path/to/data/protein.pdbqt
protein_setup: /path/to/data/grid.txt
exhaustiveness: 8
scoring: default
cnn_scoring: rescore
cnn: dense_ensemble
n_poses: 10
addH: False
ncpu: 1
seed: 0
```

**Important:** Mount directories containing input files:

- Apptainer: `-B /path/to/data:/path/to/data`
- Docker: `-v /path/to/data:/path/to/data`

## Distributed Computing

Distribute docking across multiple servers using Dask.

### Setup

**1. Create hostfile:**

With SLURM:
```bash
srun hostname | sort | uniq > NODEFILE
```

With PBS:
```bash
# Use $PBS_NODEFILE directly
```

**2. Start Dask cluster and run docking:**

```bash
dask ssh --hostfile $PBS_NODEFILE --nworkers 15 --nthreads 1 &
sleep 10
easydock -i input.smi -o output.db --program vina --config config.yml \
         --sdf --hostfile $PBS_NODEFILE --dask_report
```

Parameters:

- `--hostfile`: File with server IP addresses
- `--nworkers`: Number of workers per host (molecules docked in parallel per host)
- `--nthreads`: Can be any value; actual CPUs taken from config
- `--dask_report`: Generate HTML performance report

!!! warning "Configuration Requirements"
    - **File limit**: Increase open file limit to at least 2× total workers
    - **SSH access**: Ensure SSH connectivity with default settings between all nodes
    - **Shared filesystem**: All nodes must access the same database file

Check file limit:
```bash
ulimit -n
```

## Reusing Databases

### Create Clean Copy

Create a copy of an existing database removing all docking data:

```bash
make_clean_copy -i original.db -o clean.db
```

This preserves:
- Initialized molecules
- Stereoisomers
- Protonation states
- Setup parameters

Use the clean copy for docking with different proteins, programs or settings skipping the ligand initialization stage.

## Boron-Containing Compounds

EasyDock automatically handles boron-containing compounds for programs that don't natively support boron:

1. Replaces boron with carbon before docking
2. Returns boron atoms after docking

No special flags needed - this happens automatically as implemented in individual docking pipelines.

## Custom Configuration Arguments

Two ways to pass arguments to external programs in `config.yml`:

### Method 1: Individual Config Entries

```yaml
exhaustiveness: 8
n_poses: 5
```

!!! warning "Limited set of parameters"
    This provides access to a limited set of parameters

### Method 2: In script_file Value

```yaml
script_file: /path/to/qvina2 --exhaustiveness 8 --num_modes 5
```

!!! important
    Each argument should be specified in only ONE way, not both. Arguments in `script_file` take precedence.
