# EasyDock Documentation

[![GitHub Repo](https://img.shields.io/badge/GitHub-ci--lab--cz%2Feasydock-blue?logo=github)](https://github.com/ci-lab-cz/easydock)
[![Stars](https://img.shields.io/github/stars/ci-lab-cz/easydock?style=social)](https://github.com/ci-lab-cz/easydock/stargazers)
[![Forks](https://img.shields.io/github/forks/ci-lab-cz/easydock?style=social)](https://github.com/ci-lab-cz/easydock/network)

Welcome to EasyDock documentation! EasyDock is a fully automatic pipeline for molecular docking that provides a comprehensive solution for drug discovery workflows.

## Overview

EasyDock automates the entire docking process from ligand preparation to result analysis, supporting multiple docking programs and providing organized result storage. EasyDock does not prepare protein structures, which should be prepared for each docking program according to their requirements.

### Key Features

- **Multiple Docking Programs**: Support for Vina, Gnina/Smina, QVina, Vina-GPU and their derivatives
- **Server-Based Docking**: Containerized docking programs (CarsiDock, SurfDock, Vina-GPU) via a persistent server protocol
- **Generic Docking**: Run any external docking binary or Python script via a YAML config file, without code changes
- **Automated Preparation**: Molecule validation, salt removal, and stereoisomer enumeration
- **Flexible Protonation**: Multiple methods including MolGpKa, Uni-pKa, Chemaxon, and pkasolver
- **Container Support**: Run docking and protonation tools through Apptainer/Singularity or Docker with automatic GPU detection
- **Distributed Computing**: Scale across multiple servers using Dask
- **Database Storage**: All results organized in SQLite databases
- **Pose Quality Assessment**: PoseBusters integration (`easydock_bust`) for physics-based validation of docked poses
- **PLIF Analysis**: Protein-ligand interaction fingerprints for detailed analysis
- **Resumable Calculations**: Interrupted runs can be continued seamlessly

### Quick Start

```bash
# Create environment
conda env create -f env.yml -n easydock
# or use mamba (should be faster) 
mamba env create -f env.yml -n easydock

# Run docking
easydock -i input.smi -o output.db --program vina --config config.yml --protonation molgpka -c 4 --sdf
```

## Getting Help

- **GitHub Repository**: [ci-lab-cz/easydock](https://github.com/ci-lab-cz/easydock)
- **Issues**: Report bugs on [GitHub Issues](https://github.com/ci-lab-cz/easydock/issues)
- **Citation**: See [Citation](citation.md) page

## Documentation Structure

- **[Installation](installation.md)**: Setup instructions for all platforms
- **[Usage Guide](usage.md)**: Complete guide for running docking simulations
- **[Advanced Features](advanced.md)**: Ring sampling, distributed computing, and more
- **[Data Retrieval](data_retrieval.md)**: Extracting and analyzing results
- **[Pose Quality Assessment](posebusters.md)**: Running PoseBusters on docked poses with `easydock_bust`
- **[PLIF Analysis](plif.md)**: Protein-ligand interaction fingerprints
- **[Python API](python_api.md)**: Using EasyDock programmatically
- **[Third-party Licenses](licenses.md)**: Licenses of all integrated tools
- **[Citation](citation.md)**: Citation information
- **[Changelog](changelog.md)**: Changes in each version

### For Developers

- **[Adding a Docking Program](custom_docking.md)**: Implementing a new `xxx_dock.py` module
- **[Adding a Protonation Tool](custom_protonation.md)**: File-based, native Python, and container-based protonation conventions
- **[Server Protocol](server_protocol.md)**: Client-server docking protocol and container interface conventions
- **[Generic Docking](generic_dock.md)**: Using EasyDock with arbitrary external docking programs
