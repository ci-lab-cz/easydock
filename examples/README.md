

1. Paths in configuration files can be relative 
2. SIF-containers should be built or downloaded into a directory `containers`. If the path differs, adjust conguration files accordingly


## Prepare ligands

MolGpKa and Uni-pKa will use GPU if it will be detected.

```bash
# adjust -c to your CPU number
easydock -i 25.smi -o output.db --protonation molgpka -c 10
easydock -i 25.smi -o output.db --protonation molgpka_fix -c 10
easydock -i 25.smi -o output.db --protonation pkasolver -c 10
easydock -i 25.smi -o output.db --protonation containers/unipka.sif -c 10
```

You can make multiple copies of `output.db` and use them for docking with different programs


# Screening

Every screening can be run from scratch (ligand preparation will be executed) or starting from a previously prepared database of ligands 

## Vina (--program vina, CPU-only)

From scratch:
```bash
# adjust -c (the number of molecules docked in parallel) to your number of CPUs 
# (note that each molecule docking will run on 4 cpus, see config file)
easydock -i 25.smi -o vina.db --protonation containers/unipka.sif --program vina --config config_vina.yml -c 50
```

From `output.db`:
```bash
easydock -o output.db --program vina --config config_vina.yml -c 50
```


## Gnina (--program gnina, GPU preferable)

To run Gnina using a container provided by the Gnina developers there may be a need to specify further CUDA variables, therefore the calling command may look more complex and a user should specify explicitly all other arguments (e.g. binding paths) and use absolute paths in `config.yml`. Adjust these values to correspond to your local paths.  

From scratch:
```bash
# adjust -c (the number of molecules docked in parallel) to your number of CPUs and memory size 
# (note that each docking will load a separate gnina instance into memory)
easydock -i 25.smi -o gnina.db --protonation containers/unipka.sif --program gnina --config config_gnina.yml -c 16
```

From `output.db`:
```bash
easydock -o output.db --program gnina --config config_gnina.yml -c 16
```


## QVina (--program server, CPU-only)

```bash
# adjust -c (the number of molecules docked in parallel) to your number of CPUs 
# (note that each molecule docking will run on 4 cpus, see config file)
easydock -i 25.smi -o qvina.db --protonation containers/unipka.sif --program server --config config_qvina.yml -c 35
```

From `output.db`:
```bash
easydock -o output.db --program server --config config_qvina.yml -c 35
```


## Vina-GPU (--program server, GPU-only)

From scratch:
```bash
easydock -i 25.smi -o vinagpu.db --protonation containers/unipka.sif --program server --config config_vinagpu.yml
```

From `output.db`:
```bash
easydock -o output.db --program server --config config_vinagpu.yml
```


## CarsiDock (--program server, GPU-only)

From scratch:
```bash
# based on our experience protonation does not much affect docking with CarsiDock
easydock -i 25.smi -o carsi.db --protonation containers/unipka.sif --program server --config config_carsidock.yml
```

From `output.db`:
```bash
easydock -o output.db --program server --config config_carsidock.yml
```


## SurfDock (--program server, GPU preferable)

From scratch:
```bash
# based on our experience protonation does not much affect docking with SurfDock
easydock -i 25.smi -o surf.db --protonation containers/unipka.sif --program server --config config_surfdock.yml
```

From `output.db`:
```bash
easydock -o output.db --program server --config config_surfdock.yml
```
