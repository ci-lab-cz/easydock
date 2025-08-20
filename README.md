# EasyDock - Python module to automate molecular docking

## Table of content
- [Installation](#installation)
  - [Dependencies](#dependencies)
- [Description](#description)
  - [Features](#features)
  - [Pipeline](#pipeline)
- [Examples](#examples)
  - [Initialization of a database](#initialization-of-a-database)
  - [Docking from command line](#docking-from-command-line)
    - [Vina](#vina)
    - [Gnina / Smina](#gnina--smina)
    - [Vina-GPU / QVina-GPU / QVinaW-GPU](#vina-gpu--qvina-gpu--qvinaw-gpu)
    - [QVina2 / Qvina-W](#qvina2--qvina-w)
  - [Sampling of saturated rings](#sampling-of-saturated-rings)
  - [Docking using multiple servers](#docking-using-multiple-servers)
  - [Docking from Python](#docking-from-python)
  - [Data retrieval from the output database](#data-retrieval-from-the-output-database)
    - [Examples](#examples-1)
- [Customization](#customization)
- [Notes](#notes)
  - [Protonation notes](#protonation-notes)
  - [Multiple CPUs](#multiple-cpus)
- [Changelog](#changelog)
- [License](#licence)
- [Citation](#citation)

## Installation

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
conda install -c conda-forge python=3.11 rdkit numpy==1.26
conda install -c conda-forge scipy dask distributed scikit-learn
# for MacOS only (additional installations)
conda install -c conda-forge boost swig
```

- from pypi

```
pip install paramiko vina
# has fixes with rdkit.six import

#recommend to install the dev version, since it implements a fix making double bound non-rotatable, unlike the current release version
pip install git+https://github.com/forlilab/Meeko.git@develop
pip install gemmi

#the current release version
pip install git+https://github.com/forlilab/Meeko.git
```

- (optional) installation of `gnina/smina` is described at https://github.com/gnina/gnina
- (optional) installation of QVina2 and QVina-W is described at https://github.com/QVina/qvina
- (optional) installation of Vina-GPU, QVina2-GPU or QVinaW-GPU is described at https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1

- (recommended) installation of MolGpKa for molecule protonation
```
# pip package
pip install git+https://github.com/ci-lab-cz/MolGpKa.git
pip install torch_geometric
pip install torch==2.2 torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
# for Linux
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.2.0+cpu.html
# for MacOS (may also works for Linux)
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-2.2.0+cpu.html 
```
- (not recommended) installation of a `pkasolver` fork to enable protonation with this tool. We recommend to install and use CPU-based version, however, one may try a GPU supported one (look at the dependencies in `pkasolver` package).
```
pip install torch==1.13.1+cpu  --extra-index-url https://download.pytorch.org/whl/cpu
pip install torch-geometric==2.0.1
pip install torch_scatter==2.1.1+pt113cpu -f https://data.pyg.org/whl/torch-1.13.1%2Bcpu.html
pip install torch_sparse==0.6.17+pt113cpu -f https://data.pyg.org/whl/torch-1.13.1%2Bcpu.html
pip install torch_spline_conv==1.2.2+pt113cpu -f https://data.pyg.org/whl/torch-1.13.1%2Bcpu.html
pip install cairosvg svgutils
pip install git+https://github.com/Feriolet/dimorphite_dl.git
pip install git+https://github.com/DrrDom/pkasolver.git
```

## Description

The fully automatic pipeline for molecular docking.  

An important feature, once the database is initialized it will store all command line arguments and input files as a part of the setup table. This will allow to easy rerun interrupted calculations and keep all data in one place for consistency, but this also makes impossible to change some settings afterwards by providing other values in command line arguments. The existing database will never be overwritten. Therefore, if the database was initialized wrongly, it should be deleted before rerun the command line script. Alternatively, incorrect values of arguments can be edited directly in the database.  

### Features

- the major script `easydock` supports docking with `vina`, `gnina` (`gnina` also supports `smina` and its custom scoring functions), `qvina` (`qvina` is a collective term to use either QVina2 or QVina-W), `vina-gpu` (`vina-gpu` is a collective term to use any of Vina-GPU, QVina2-GPU or QVinaW-GPU)
- can be used as a command line utility or imported as a python module
- if input molecules are 3D, these conformations will be used as starting ones for docking (enable usage of external conformer generators)
- input molecules are checked for salts and attempted to fix by SaltRemover 
- stereoisomers can be enumerated for unspecified chiral centers and double bonds (since some compounds may require very long runtimes, the maximum runtime for individual molecules was set to 300 sec)
- several protonation options: molgpka, chemaxon and pkasolver (check notes below). If omitted input protonation state will be used (enables usage of external protonation tools)
- docking of compounds with saturated rings can be enhanced by additional sampling of starting ring conformers and only the best one is stored to the database
- supports distributed computing using `dask` library
- supports docking of boron-containing compounds using `vina` and `smina` (boron is replaced with carbon before docking and returned back)
- all outputs are stored in an SQLite database
- interrupted calculations can be continued by invoking the same command or by supplying just a single argument (`--output`) - the existing output database
- `get_sdf_from_dock_db` is used to extract data from output DB
- all command line arguments and input files are stored in the setup table and the majority of those parameters cannot be changed later. This will prevent losing input settings. If some changes should be made after DB was created a direct editing the DB can be a solution   

### Pipeline

The pipeline consists of two major parts which can be run separately or simultaneously:
1. Initialization of database, which includes:
- input SMILES are converted in 3D by RDKit, if input is 3D structures in SDF their conformations wil be taken as starting without changes.
- compounds having salts are stripped, if this fails the whole compound will be omitted for docking reporting to STDERR 
- up to a specified number of stereoisomers are enumerated for molecules with undefined chiral centers or double bond configurations (by default 1 random but reproducible stereoisomer is generated)
- ligands are protonated by MolGpKa/Chemaxon/pKasolver at pH 7.4 and the most stable tautomers are generated (optional, requires a Chemaxon license)
2. Docking step includes:
- molecules are converted in PDBQT format using Meeko
- docking with `vina`/`gnina`/`qvina`/`vina-gpu`
- top docked poses are converted in MOL format and stored into output DB along with docking scores

### Notes

- These two parts of the pipeline allows to create a DB and reuse it for docking with different proteins/settings/etc.  
- There is also a script `make_clean_copy` which creates a copy of an existing DB removing all docking data to use it for docking with different proteins/settings/etc.        
- Protonation with MolGpKa will run in a single cpu mode, but it will use 25-50% of available cores.   

## Examples

### Initialization of a database

This will create a DB with checked molecules using 4 cores. If `--protonation` argument was not used molecules will keep their input protonations states (this enables docking of molecules protonated by an external tool)
```
easydock -i input.smi -o output.db -c 4
```
Initialize DB and protonate molecules with MolGpKa
```
easydock -i input.smi -o output.db -c 4 --protonation molgpka
```

### Docking from command line

#### Vina

To run a both steps of the full pipeline: initialization and docking.  
Docking using `vina` takes input SMILES and a config file. Ligands will not be protonated. 4 CPU cores will be used (4 molecules will be docked in parallel). When docking will finish an SDF file will be created with top docking poses for each ligand. 
```
easydock -i input.smi -o output.db --program vina --config config.yml -c 4 --sdf
``` 

The same command for a previously initialized DB. Supplying arguments relevant for the initialization stage will have no effect after DB was created
```
easydock -o output.db --program vina --config config.yml -c 4 --sdf
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

NOTE: ncpu argument in `easydock` and `config.yml` has different meaning. In `easydock` it means the number of molecules docked in parallel. In `config.yml` it means the number of CPUs used for docking of a single molecule. The product of these two values should be equal or a little bit more than the number of CPUs on a computer.

#### Gnina / Smina

The same but using `gnina`
```
easydock -i input.smi -o output.db --program gnina --config config.yml -c 4 --sdf
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

To use `gnina` (or other tools) from within a `docker/apptainer` container it is necessary to change the `script_file` argument. There is also a need to bind a local directory to a directory inside a container using the same path (in `docker` there is `-v` argument for that).
```
script_file: apptainer exec -B /path/to/dir:/path/to/dir container.sif gnina
protein: /path/to/dir/protein.pdbqt
protein_setup: /path/to/dir/grid.txt
exhaustiveness: 8
scoring: default
cnn_scoring: rescore
cnn: dense_ensemble
n_poses: 10
addH: False
ncpu: 1
seed: 0
```

#### Vina-GPU / QVina-GPU / QVinaW-GPU

To use any of Vina-GPU, QVina2-GPU or QVinaW-GPU use `--program vina-gpu`:

```
easydock -i input.smi -o output.db --program vina-gpu --config config.yml--sdf
```

Example of config.yml for `vina-gpu` docking:  
```
script_file: /path/to/bin/dir/AutoDock-Vina-GPU-2-1 --opencl_binary_path /path/to/bin/dir/
protein: /path/to/dir/protein.pdbqt
protein_setup: /path/to/dir/grid.txt
n_poses: 3
```
`--opencl_binary_path` is necessary to get access to OpenCL macros.
`thread` argument is set to 8000 by default, this can be changed by explicitly passing its value in `config.yml`.

Choosing other programs `QuickVina2-GPU-2-1` or `QuickVina-W-GPU-2-1` in `config.yml` will enable docking with them. All other arguments are the same.

#### QVina2 / QVina-W

To use QVina2 or QVina-W use `--program qvina`

```
easydock -i input.smi -o output.db --program qvina --config config.yml--sdf
```

Example of config.yml for `qvina` docking:  
```
script_file: /path/to/bin/dir/qvina2
protein: /path/to/dir/protein.pdbqt
protein_setup: /path/to/dir/grid.txt
exhaustiveness: 8
n_poses: 3
ncpu: 2
seed: 0
```

Specify path to either QVina2 or QVina-W to use a particular program.

#### Config.yml

There are two ways to pass arguments to programs which are executed by external binaries (currently these are all programs except vina). The first way is to specify arguments and their values as individual entries in `config.yml`. However, not all arguments can be set in this way. The second approach is specify constant arguments in the `script_file` value. For example, if you want to customize `--energy_range` for `qvina` use `script_file: /path/to/bin/dir/qvina2 --energy_range 5`. Each argument should be specified in either way, not both.

### Sampling of saturated rings

Since rings are considered rigid, to improve docking of compounds containing saturated rings a set of starting conformers can be sampled and used for docking (enabled by `--ring_sample` argument). Only the conformer resulted int he best docking score will be stored in the database. This will increase computational complexity proportionally and may result in much longer simulation times especially for compounds containing multiple saturated rings.  
```
easydock -i input.smi -o output.db --program vina --config config.yml -c 4 --ring_sample 
```

### Docking using multiple servers

To distribute docking over multiple servers one have to start dask cluster and call the script as follows. 

```bash
dask ssh --hostfile $PBS_NODEFILE --nworkers 15 --nthreads 1 &
sleep 10
easydock -i input.smi -o output.db --program vina --config config.yml --sdf --hostfile $PBS_NODEFILE --dask_report
```
`$PBS_NODEFILE` is a file containing list of IP addresses of servers. The first one from the list will be used by a dask scheduler, but it will also participate in computations.  
To create this file with SLURM one may use the following command `srun hostname | sort | uniq > NODEFILE` 

`--nworkers` is the number of workers per host. This is the number of molecules which are docked in parallel on a single host.

`--nthreads` can be any value. The number of CPUs used for docking of a single molecule will be taken from `config.yml`.
  
`--dask_report` argument will create at the end of calculations a html-file with performance report (may be useful to tweak docking parameters).  
  
**Important setup issues**:
- the limit of open files on every server should be increased to the level at least twice the total number of requested workers (file streams are used for internode communication by dask).
- all nodes should be accessible by SSH with default settings 

### Docking from Python

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

### Data retrieval from the output database

1. Using `--sdf` argument of the main script `easydock` will return top poses with docking scores. If there were several enumerated stereoisomers, it will return the pose and the score of the best scoring stereoisomer only.  

2. Using `get_sdf_from_dock_db` script. it has a rich set of settings and can return SDF as well as SMILES files. The only restriction it cannot currently return the best pose among enumerated stereoisomers. In this case it is advised to use the previous option and invoke `easydock -o database.db --sdf` on the database with docked molecules.

#### Examples

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

## Customization

To implement support of a custom docking program one should implement a function like `mol_dock` which will take as input an RDKit mol object (named molecule) and a yml-file with all docking parameters. The function should run a command line script/utility and return back a tuple of a molecule name and a dictionary of parameters and their values which should be stored in DB (parameter names should be exactly the same as corresponding field names in DB). For examples, please look at `mol_dock` functions in `vina_dock` or `gnina_dock`.

Details on implementation of support of other protonation tool are given in the `protonation.py` module. 

## Notes

### Protonation notes

There are two integrated open-source approaches (molgpka and pkasolver) and one commercial (chemaxon).
1. Chemaxon is pretty reliable, however it requires a paid license.
2. MolGpKa is a model trained on predictions of Chemaxon and thus aligns well with it. Protonation states of each atom are chosen based on the predicted pKa and pKb values. However, there are certain issues. Some issues were fixed by post-processing SMARTS patterns (avoid protonation of amide groups, etc). However, some issues are difficult to fix (e.g. missing protonation of aliphatic amines, piperizines, etc).  
3. pkasolver enumerates protonation states and the form closest to pH 7.4 is selected as a relevant one. In some cases it may return invalid SMILES, e.g. `O=C(N1CCN(CC1)C(=O)C=2C=CC=C(C#CC3CC3)C2)C=4NN=C5CCCC45 -> O=C(c1cccc(C#CC2CC2)c1)N1CC[NH](C(=O)c2[nH]nc3c2CCC3)CC1`, which will be skipped and a corresponding warning message will appear. It has many issues with protonated forms (very frequent protonation of amide groups, etc). Therefore, we do not recommend its usage.

### Multiple CPUs

Please pay attention for `--ncpu` argument if you use `--protonation pkasolver`. For `ncpu` > 1 it may result in some errors. Please report this issue. 

## Changelog

[changelog.md](changelog.md)

## Licence
BSD-3

## Citation
Minibaeva, G.; Ivanova, A.; Polishchuk, P.,  
EasyDock: customizable and scalable docking tool.  
*Journal of Cheminformatics* **2023**, 15 (1), 102.  
https://doi.org/10.1186/s13321-023-00772-2
