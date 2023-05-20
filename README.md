# EasyDock - Python module to automate molecular docking

### Installation

```
pip install easydock
```
or the latest version from github
```
pip install git+https://github.com/ci-lab-cz/easydock.git
```

### Dependencies

from conda
```
conda install -c conda-forge python=3.9 numpy=1.20 rdkit scipy dask distributed
```

from pypi
```
pip install vina meeko
```

Installation of gnina is described at https://github.com/gnina/gnina

### Description

Fully automatic pipeline for molecular docking.  

Features:
- the major script `run_dock` supports docking with `vina` and `gnina` (`gnina` also supports `smina` and its custom scoring functions)
- can be used as a command line utility or imported as a python module
- supports distributed computing using `dask` library
- supports docking of boron-containing compounds using `vina` and `smina` (boron is replaced with carbon before docking and returned back)
- all outputs are stored in an SQLite database
- interrupted calculations can be restarted by invoking the same command or by supplying just a single argument - the existing output database
- `get_sdf_from_dock_db` is used to extract data from output DB

Pipeline:
- input SMILES are converted in 3D by RDKit, if input is 3D structures in SDF their conformations wil be taken as starting without changes.
- ligands are protonated by chemaxon at pH 7.4 and the most stable tautomers are generated (optional, requires a Chemaxon license)
- molecules are converted in PDBQT format using Meeko
- docking with `vina`/`gnina`
- output poses are converted in MOL format and stored into output DB along with docking scores

### Example

##### Docking from command line

Docking using `vina` takes input SMILES and a config file. Ligands will not be protonated with Chemaxon, so their supplied charged states will be used. 4 CPU cores will be used. When docking will finish an SDF file will be created with top docking poses for each ligand. 
```
run_dock -i input.smi -o output.db --program vina -c config.yml --no_protonation -c 4 --sdf
``` 

Example of config.yml for `vina` docking  
```
protein: /path/to/protein.pdbqt
protein_setup: /path/to/grid.txt
exhaustiveness: 8
seed: 0
n_poses: 5
ncpu: 5
```

The same but using `gnina`
```
run_dock -i input.smi -o output.db --program gnina -c config.yml --no_protonation -c 4 --sdf
``` 

Example of config.yml for `gnina` docking  
```
script_file: /path/to/gnina_executable
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

To use `smina` invoke `gnina` as shown above and make corresponding changes in config.yml
```
script_file: /path/to/gnina_executable
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

##### Docking using multiple servers

To distribute docking over multiple servers one have to start dask cluster and call the script

```bash
dask ssh --hostfile $PBS_NODEFILE --nworkers 15 --nthreads 1 &
sleep 10
run_dock -i input.smi -o output.db --program vina -c config.yml --no_protonation -c 4 --sdf --hostfile $PBS_NODEFILE --dask_report
```
`$PBS_NODEFILE` is a file containing list of IP addresses of servers. The first one from the list will be used by a dask scheduler, but it will also participate in computations.
  
`--dask_report` argument will create at the end of calculations an html-file with performance report (may be useful to tweak docking parameters).
  
**Important setup issue** - the limit of open files on every server should be increased to the level at least twice the total number of requested workers (file streams are used for inter-node communication by dask).

##### Data retrieval from the output database

To extract data from the database one may use the script `get_sdf_from_dock_db`.

Extract top poses with their scores (additional information in DB fields can be extracted only for the top poses):
```
get_sdf_from_dock_db -i output.db -o output.sdf --fields docking_score
```
Retrieve second poses for compounds `mol_1` and `mol_4` in SDF format:
```
get_sdf_from_dock_db -i output.db -o output.sdf -d mol_1 mol_4 --poses 2 
```
Instead of a list of ids a text file can be supplied as an argument `-d`.

Retrieve top poses for compounds with docking score less then -10:
```
get_sdf_from_dock_db -i output.db -o output.sdf --fields docking_score --add_sql 'docking_score < -10' 
```

##### Docking from Python

Dock a list of molecules on a local computer. Import `mol_dock` function from a corresponding submodule.
```python
from easydock.run_dock import docking
from easydock.vina_dock import mol_dock
# from easydock.gnina_dock import mol_dock  # <- enable gnina docking
from rdkit import Chem

smiles = ['CC(=O)O', 'NCC(=O)O', 'NC(C)C(=O)O']
mols = [Chem.MolFromSmiles(smi) for smi in smiles]

# assign names, because this will be an identifier of docking outputs of a molecule 
for mol, smi in zip(mols, smiles):
    mol.SetProp('_Name', smi)

for mol_id, res in docking(mols, dock_func=mol_dock, dock_config='config.yml', ncpu=4):
    print(mol_id, res)
```

##### Customization

To implement support of a custom docking program one should implement a function like `mol_dock` which will take as input an RDKit mol object (named molecule) and an yml-file with all docking parameters. The function should run a command line script/utility and return back a tuple of a molecule name and a dictionary of parameters and their values which should be stored in DB (parameter names should be exactly the same as corresponding field names in DB). For examples, please look at `mol_dock` functions in `vina_dock` or `gnina_dock`.

### Changelog

**0.2.0**
- the stable version with multiple fixes and updates
- dask library was fully integrated and tested
- API was redesigned
- docking of boron-containing compounds was implemented (Vina, smina)
- a function to predict docking runtime was introduced   

**0.1.2**
- (bugfix) docking of macrocycles is rigid (in future may be changed)

### Licence
BSD-3
