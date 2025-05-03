**1.0.0**
- separate initialization of a database (ligand preparation) and docking steps
- create a script to clear docking data from a database to reuse it for other docking (different protein, settings, etc)
- speed up initialization of a database
- integration of MolGpKa for molecule protonation and implement post-processing SMARTS patterns to fix errors 
- implement optional sampling of saturated ring systems to improve docking results 
- fix some errors in pkasolver by applying post-processing SMARTS patterns
- make clear command argument description and grouping

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
