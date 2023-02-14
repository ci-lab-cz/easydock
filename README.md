# Python scripts to automate molecular docking

### Installation

```
pip install moldock
```
or the latest version
```
pip install git+https://github.com/ci-lab-cz/docking-scripts.git
```

### Dependencies

from conda
```
conda install -c conda-forge python=3.9 numpy=1.20 rdkit scipy dask distributed vina
```

from pypi
```
pip install meeko
```

Installation of gnina is described at https://github.com/gnina/gnina

### Description

Fully automatic pipeline for molecular docking.
- two major scripts `vina_dock` and `gnina_dock` which support docking with `vina` and `gnina` (`gnina` also supports `smina` and its custom scoring functions)
- can be used as command line scripts or imported as a python module
- support distributed computing using `dask` library
- `get_sdf_from_db` is used to extract data from output DB 

Pipeline:
- input SMILES are converted in 3D by RDKit embedding, if input is 3D structures in SDF their conformations wil be taken as starting without changes.
- ligands are protonated by chemaxon at pH 7.4 and the most stable tautomers are generated (optional, requires a Chemaxon license)
- molecules are converted in PDBQT format
- docking with `vina`/`gnina`
- output poses are converted in MOL format and stored into output DB along with docking scores

### Example

Both scripts `vina_dock` and `gnina_dock` have similar common arguments.

Docking using input SMILES, prepared protein and config files. Ligands will not be protonated with Chemaxon, so their supplied charged states will be used. 4 CPU cores will be used. When docking will finish an SDF will be created with top docking poses for each ligand. 
```
vina_dock -i input.smi -o output.db -p protein.pdbqt -s vina_config --no_protonation -c 4 --sdf 
``` 

Retrieve second poses for compounds `mol_id_1` and `mol_id_2` with their docking scores in SDF format:
```
get_sdf_from_db -i output.db -o out.sdf -d mol_id_1,mol_id_4 --fields docking_score --poses 2 
```
Instead of a comma-separated list of ids a text file can be supplied as an argument `-d`.

Retrieve top poses for compounds with docking score less then -10:
```
get_sdf_from_db -i output.db -o out.sdf --fields docking_score --poses 1 --add_sql 'docking_score < -10' 
```

### Changelog

**0.1.2**
- (bugfix) docking of macrocycles is rigid (in future may be changed)

### Licence
CC BY-NC-SA 4.0
