# EasyDock Documentation

Welcome to EasyDock documentation! EasyDock is a fully automatic pipeline for molecular docking that provides a comprehensive solution for drug discovery workflows.

## Overview

EasyDock automates the entire docking process from molecule preparation to result analysis, supporting multiple docking programs and providing organized result storage.

### Key Features

- **Multiple Docking Programs**: Support for Vina, Gnina/Smina, QVina, Vina-GPU and their derivatives
- **Automated Preparation**: Molecule validation, salt removal, and stereoisomer enumeration
- **Flexible Protonation**: Multiple methods including MolGpKa, Uni-pKa, Chemaxon, and pkasolver
- **Distributed Computing**: Scale across multiple servers using Dask
- **Database Storage**: All results organized in SQLite databases
- **PLIF Analysis**: Protein-ligand interaction fingerprints for detailed analysis
- **Resumable Calculations**: Interrupted runs can be continued seamlessly

### Quick Start

```bash
# Create environment
conda env create -f env.yml

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
- **[PLIF Analysis](plif.md)**: Protein-ligand interaction fingerprints
- **[Configuration](configuration.md)**: Parameter reference
- **[Python API](python_api.md)**: Using EasyDock programmatically
- **[Customization](customization.md)**: Extending EasyDock functionality
- **[Troubleshooting](troubleshooting.md)**: Solutions to common problems
