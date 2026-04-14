# Installation

## Quick Installation

### Using Conda Environment (Recommended)

The easiest way to install EasyDock with all dependencies:

```bash
conda env create -f env.yml
```

This installs EasyDock with Vina (docking) and MolGpKa (protonation).

## Manual Installation

### Install EasyDock Package

From PyPI:
```bash
pip install easydock
```

From GitHub (latest development version):
```bash
pip install git+https://github.com/ci-lab-cz/easydock.git
```

### Required Dependencies

#### From Conda

```bash
conda env create -n easydock
conda install -c conda-forge python=3.11 rdkit numpy==1.26
conda install -c conda-forge scipy dask distributed scikit-learn
```

Additional for macOS:
```bash
conda install -c conda-forge boost swig
```

#### From PyPI

```bash
pip install paramiko vina prolif gemmi
```

Meeko (recommended development version):
```bash
pip install git+https://github.com/forlilab/Meeko.git@develop
```

The development version includes fixes that make double bonds non-rotatable, improving accuracy.

Or current release version:
```bash
pip install git+https://github.com/forlilab/Meeko.git
```

## Optional Docking Programs

### Gnina/Smina

Installation instructions at [https://github.com/gnina/gnina](https://github.com/gnina/gnina)

### QVina2 and QVina-W

Installation instructions at [https://github.com/QVina/qvina](https://github.com/QVina/qvina)

### Vina-GPU Family

Vina-GPU, QVina2-GPU, and QVinaW-GPU installation at [https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1)

## Optional Protonation Tools

### MolGpKa (Recommended)

```bash
pip install git+https://github.com/ci-lab-cz/MolGpKa.git
pip install torch_geometric
pip install torch==2.2 torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
```

**For Linux:**
```bash
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv \
    -f https://data.pyg.org/whl/torch-2.2.0+cpu.html
```

**For macOS:**
```bash
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv \
    -f https://data.pyg.org/whl/torch-2.2.0+cpu.html
```

### pkasolver

```bash
pip install torch==1.13.1+cpu --extra-index-url https://download.pytorch.org/whl/cpu
pip install torch-geometric==2.0.1
pip install torch_scatter==2.1.1 -f https://data.pyg.org/whl/torch-1.13.1%2Bcpu.html
pip install torch_sparse==0.6.17 -f https://data.pyg.org/whl/torch-1.13.1%2Bcpu.html
pip install torch_spline_conv==1.2.2 -f https://data.pyg.org/whl/torch-1.13.1%2Bcpu.html
pip install cairosvg svgutils
pip install git+https://github.com/Feriolet/dimorphite_dl.git
pip install git+https://github.com/DrrDom/pkasolver.git
```

!!! warning "pkasolver Issues"
    pkasolver is not recommended due to frequent errors and invalid SMILES generation. Use MolGpKa or Uni-pKa instead.

## Container Support

EasyDock supports running external tools through Singularity/Apptainer containers.

### Linux

Install Singularity or Apptainer following their official documentation.

### macOS

EasyDock supports Docker natively on macOS. Install Docker and ensure CLI access without sudo — no additional setup is required.

#### Running SIF Containers via Apptainer Proxy (Advanced)

If you need to run Apptainer SIF containers on macOS, you can set up a Docker-based Apptainer proxy. Create a file named `Dockerfile`:

```dockerfile
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    wget build-essential squashfs-tools uidmap fuse3 git \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/apptainer/apptainer/releases/download/v1.3.4/apptainer_1.3.4_amd64.deb \
    && dpkg -i apptainer_1.3.4_amd64.deb

ENTRYPOINT ["apptainer"]
```

Switch to the directory with `Dockerfile` and build the container:

```bash
docker build -t apptainer:latest --platform=linux/amd64 .
```

EasyDock will automatically use the `apptainer:latest` image as a proxy to run SIF containers.

### Windows

Install EasyDock in Windows Subsystem for Linux (WSL) to use container features.

### Pre-built Containers

| Container | Download | Description |
|-----------|----------|-------------|
| unipka.sif | [Zenodo](https://doi.org/10.5281/zenodo.17506577) | Uni-pKa model for protonation |

**Usage Example:**
```bash
easydock -i input.smi -o output.db -c 4 --protonation /path/to/unipka.sif
```

### Building Containers from Source

Definition files and Dockerfiles for all EasyDock containers are in the `containers/` directory of the repository. Apptainer (`.def`) and Docker images can be built interchangeably — both expose the same command interface.

#### Uni-pKa (protonation)

```bash
cd containers/unipka

# Apptainer
apptainer build unipka.sif unipka.def

# Docker
docker build -t unipka .
```

#### CarsiDock (server-based docking)

The build downloads model checkpoints (~1 GB) from Google Drive and GitHub. An internet connection is required.

```bash
cd containers/carsidock

# Apptainer
apptainer build carsidock.sif carsidock_server.def

# Docker
docker build -t carsidock .
```

#### Vina-GPU Server (server-based docking)

The build downloads CPU and GPU binaries and Boost from source. Requires an internet connection.

```bash
cd containers/vinagpu

# Apptainer only (no Dockerfile provided)
apptainer build vinagpu.sif vinagpu_server.def
```

## Optional Analysis Tools

### PoseBusters

PoseBusters provides physics-based validation of docked poses. Required for the `easydock_bust` tool and the `--bust` flag of `get_sdf_from_easydock`.

```bash
pip install posebusters
```

## Verification

Verify installation:

```bash
easydock --help
```

You should see the help message with available options.
