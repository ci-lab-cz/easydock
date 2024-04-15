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

- from conda
```
conda install -c conda-forge python=3.9 numpy=1.20 rdkit scipy dask distributed
```

- from pypi

```
pip install paramiko meeko vina
```

- (optional) installation of `gnina/smina` is described at https://github.com/gnina/gnina

- (optional) installation of a `pkasolver` fork to enable protonation with this tool. We recommend to install and use CPU-based version, however, one may try a GPU supported one (look at the dependencies in `pkasolver` package).
```
pip install torch==1.13.1+cpu  --extra-index-url https://download.pytorch.org/whl/cpu
pip install torch-geometric==2.0.1
pip install torch_scatter==2.1.1+pt113cpu -f https://data.pyg.org/whl/torch-1.13.1%2Bcpu.html
pip install torch_sparse==0.6.17+pt113cpu -f https://data.pyg.org/whl/torch-1.13.1%2Bcpu.html
pip install torch_spline_conv==1.2.2+pt113cpu -f https://data.pyg.org/whl/torch-1.13.1%2Bcpu.html```
pip install git+https://github.com/Feriolet/dimorphite_dl.git
pip install git+https://github.com/DrrDom/pkasolver.git
```

### Description

The fully automatic pipeline for molecular docking.  

Features:
- the major script `run_dock` supports docking with `vina` and `gnina` (`gnina` also supports `smina` and its custom scoring functions)
- can be used as a command line utility or imported as a python module
- stereoisomers can be enumerated for unspecified chiral centers and double bonds
- several protonation options: chemaxon and pkasolver (check notes below)
- supports distributed computing using `dask` library
- supports docking of boron-containing compounds using `vina` and `smina` (boron is replaced with carbon before docking and returned back)
- all outputs are stored in an SQLite database
- interrupted calculations can be continued by invoking the same command or by supplying just a single argument - the existing output database
- `get_sdf_from_dock_db` is used to extract data from output DB

Pipeline:
- input SMILES are converted in 3D by RDKit, if input is 3D structures in SDF their conformations wil be taken as starting without changes.
- up to a specified number of stereoisomers are enumerated for molecules with undefined chiral centers or double bond configurations (by default 1 random but reproducible stereoisomer is generated)
- ligands are protonated by Chemaxon/pKasolver at pH 7.4 and the most stable tautomers are generated (optional, requires a Chemaxon license)
- molecules are converted in PDBQT format using Meeko
- docking with `vina`/`gnina`
- output poses are converted in MOL format and stored into output DB along with docking scores

### Example

##### Docking from command line

Docking using `vina` takes input SMILES and a config file. Ligands will not be protonated with Chemaxon, so their supplied charged states will be used. 4 CPU cores will be used (4 molecules will be docked in parallel). When docking will finish an SDF file will be created with top docking poses for each ligand. 
```
run_dock -i input.smi -o output.db --program vina --config config.yml --no_protonation -c 4 --sdf
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

NOTE: ncpu argument in `run_dock` and `config.yml` has different meaning. In `run_dock` it means the number of molecules docked in parallel. In `config.yml` it means the number of CPUs used for docking of a single molecule. The product of these two values should be equal or a little bit more than the number of CPUs on a computer.

The same but using `gnina`
```
run_dock -i input.smi -o output.db --program gnina --config config.yml --no_protonation -c 4 --sdf
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
run_dock -i input.smi -o output.db --program vina --config config.yml --no_protonation --sdf --hostfile $PBS_NODEFILE --dask_report
```
`$PBS_NODEFILE` is a file containing list of IP addresses of servers. The first one from the list will be used by a dask scheduler, but it will also participate in computations.

`--nworkers` is the number of workers per host. This is the number of molecules which are docked in parallel on a single host.

`--nthreads` can be any value. The number of CPUs used for docking of a single molecule will be taken from `config.yml`.
  
`--dask_report` argument will create at the end of calculations a html-file with performance report (may be useful to tweak docking parameters).  
  
**Important setup issue** - the limit of open files on every server should be increased to the level at least twice the total number of requested workers (file streams are used for internode communication by dask).

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

Retrieve top poses for compounds with docking score less than -10:
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

##### Retrieval output poses

1. Using `--sdf` option of the main script `run_dock` will return top poses with docking scores. If there were several enumerated stereoisomers, it will return the pose and the score of the best scoring stereoisomer only.
2. Using `get_sdf_from_dock_db` script. it has a rich set of settings and can return SDF as well as SMILES files. The only restriction it cannot currently return the best pose among enumerated stereoisomers. In this case it is advised to use the previous option and invoke `run_dock -o database.db --sdf` on the database with docked molecules.

##### Customization

To implement support of a custom docking program one should implement a function like `mol_dock` which will take as input an RDKit mol object (named molecule) and a yml-file with all docking parameters. The function should run a command line script/utility and return back a tuple of a molecule name and a dictionary of parameters and their values which should be stored in DB (parameter names should be exactly the same as corresponding field names in DB). For examples, please look at `mol_dock` functions in `vina_dock` or `gnina_dock`.

##### Protonation notes

pkasolver enumerated protonation states and the closest to pH 7.4 is chosen. In some cases it may return invalid SMILES, e.g. `O=C(N1CCN(CC1)C(=O)C=2C=CC=C(C#CC3CC3)C2)C=4NN=C5CCCC45 -> O=C(c1cccc(C#CC2CC2)c1)N1CC[NH](C(=O)c2[nH]nc3c2CCC3)CC1`, which will be skipped and a corresponding warning message will appear.

Please note, that protonation states generated with `pkasolver` were not validated. So, check protonation states.

### Changelog

**0.3.0**
- add optional enumeration of stereoisomers. This partially breaks compatibility - docking of molecules in databases which were created by the previous version cannot be continued with this version. Everything else including API is compatible [#21](https://github.com/ci-lab-cz/easydock/pull/21)
- fix minor errors in retrieving non-top poses by `get_sdf_from_dock_db`

**0.2.9**
- fix extraction of docking scores for gnina outputs (critial fix) [#23](https://github.com/ci-lab-cz/easydock/pull/23)
- fix database update freezing upon errors occurred with docking individual docking programs [#22](https://github.com/ci-lab-cz/easydock/issues/22)
- implement to keep the order of output molecules in get_sdf_from_dock_db script if an argument -d is used

**0.2.8**
- conversion of PDBQT to Mol by means of Meeko (improvement) [#19](https://github.com/ci-lab-cz/easydock/pull/19)
- clarify installation instructions

**0.2.7**
- add an optional UNIQUE constraint on SMILES field in the main table on database creation (currently duplicates are not removed)

**0.2.6**
- fix compatibility issue with meeko version 0.5.0

**0.2.5**
- fix input argument type
- update examples and citation

**0.2.4**
- close pool explicitly to solve issue with multiprocessing
- replace subprocess calls with run
- explicitly set types of command line arguments which a file paths (solve issue with relative paths)

**0.2.3**
- improve descriptions of examples on README
- catch all exceptions in conversion of PDBQT to Mol
- move DB related functions to a new database.py module
- use SMILES temporary file to protonate molecules with cxcalc
- add functions to get molecules from DB in Python (get_mols, select_from_db)

**0.2.2**
- fix bug with continuation of calculations after db was transferred to other machine
- restrict precedence of command line arguments over arguments restored from DB only to specific ones (output, hostfile, dask_report, ncpu, verbose)

**0.2.1**
- fix treatment of molecule ids in get_sdf_from_dock_db
- change installation instructions, vina must be installed from sources
- add argument no_tautomerization to disable tautomerization during protonation
- (critical) fix conversion of PDBQT to Mol which could not assign bond orders and returned molecules with only single bonds 

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

### Citation
Minibaeva, G.; Ivanova, A.; Polishchuk, P.,  
EasyDock: customizable and scalable docking tool.  
*Journal of Cheminformatics* **2023**, 15 (1), 102.  
https://doi.org/10.1186/s13321-023-00772-2
