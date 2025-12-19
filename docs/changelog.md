**1.2.0**
- implement support of singularity/apptainer containers incorporating protonation tools on Linux and MacOS
- integrate uni-pka protonation as a container
- integrate protein-ligand fingerprint detection for docking poses of ligands based on ProLIF
- add env.yml to make setup and installation easier
- the console script get_sdf_from_dock_db was renamed to get_sdf_from_easydock
- documentation was created

**1.1.0**
- integrate QVina2, QVina-W, Vina-GPU, QVina2-GPU, QVinaW-GPU
- make easydock running on Windows, however, this will work only with docking programs which can be run on Windows (e.g. Python Vina cannot be installed on Windows, but if you have a binary you can run docking with it)
- add the lowest energy conformer to a set of representative structures selected during the ring sampling
- improve SMARTS fix rules for MolGpKa protonation prediction
- add a header to output of get_sdf_from_dock_db if output is SMILES format
- fix docking timestamp (now the local time is used)
- update README with more examples (e.g. use docking programs from containers) and improved the description
- improve logging

**1.0.0**
- separate initialization of a database (ligand preparation) and docking steps [#43](https://github.com/ci-lab-cz/easydock/pull/43)
- create a script to clear docking data from a database to reuse it for other docking (different protein, settings, etc)
- speed up initialization of a database
- integration of MolGpKa for molecule protonation and implement post-processing SMARTS patterns to fix errors [#49](https://github.com/ci-lab-cz/easydock/pull/49) 
- implement optional sampling of saturated ring systems to improve docking results [#34](https://github.com/ci-lab-cz/easydock/pull/34) 
- fix some errors in pkasolver by applying post-processing SMARTS patterns
- introduced logging, if a log-file is not specified the output will be redirected to stdout 
- improved command argument descriptions and grouping arguments

**0.3.2**
- catch and report all exceptions when reconstruct 3D geometry of protonated molecules
- support of temporary directories [#26](https://github.com/ci-lab-cz/easydock/issues/26) [#36](https://github.com/ci-lab-cz/easydock/pull/36)
- speed up pkasolver [#39](https://github.com/ci-lab-cz/easydock/pull/39) [#35](https://github.com/ci-lab-cz/easydock/issues/35)
- stereo_id column default value set to 0
- fix reading setup table if some fields are empty

**0.3.1**
- add pkasolver as a protonation tool ([@Feriolet](https://github.com/Feriolet)) [#17](https://github.com/ci-lab-cz/easydock/issues/17)
- add preparation step which tries to strip salts [#30](https://github.com/ci-lab-cz/easydock/issues/30)
- database structure was modified to store input molecules before pre-preprocessing [#31](https://github.com/ci-lab-cz/easydock/pull/31)
- speed of database initialization was improved [#29](https://github.com/ci-lab-cz/easydock/pull/29)

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
