<p align="left">
  <img src="docs/assets/light.svg" alt="EasyDock logo" width="300">
</p>

# EasyDock - Python module to automate molecular docking

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
conda env create -f env.yml -n easydock
# or use mamba (should be faster) 
mamba env create -f env.yml -n easydock

# Run docking
easydock -i input.smi -o output.db --program vina --config config.yml --protonation molgpka -c 4 --sdf
```

## Documentation

<https://easydock.readthedocs.io/en/latest/>

## Licence
BSD-3

## Citation
Minibaeva, G.; Ivanova, A.; Polishchuk, P.,  
EasyDock: customizable and scalable docking tool.  
*Journal of Cheminformatics* **2023**, 15 (1), 102.  
https://doi.org/10.1186/s13321-023-00772-2
