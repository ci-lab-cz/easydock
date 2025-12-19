# Protein-Ligand Interaction Fingerprints

EasyDock can calculate protein-ligand interaction fingerprints (PLIFs) using ProLIF and compute similarity to reference interactions.

## Calculate PLIFs

### From Database

Calculate PLIFs for docked molecules:

```bash
easydock_plif -i output.db -p protein.pdb -c 4
```

**Parameters:**

- `-i`: Database with docked molecules
- `-p`: Protein PDB file with all hydrogens
- `-c`: Number of CPU cores

**Output:** PLIFs are stored in the database for later retrieval.

### Specific Poses

By default, only top poses are processed. Specify other poses:

```bash
easydock_plif -i output.db -p protein.pdb --poses 1 2 3 -c 4
```

### From SDF Files

Calculate PLIFs directly from SDF (useful for reference ligands):

```bash
easydock_plif -i reference.sdf -p protein.pdb -o plif_output.txt -c 4
```

## Retrieve PLIFs

Extract calculated PLIFs from database:

```bash
easydock_plif -i output.db -o plif_output.txt
```

**Output format:** Tab-separated file with PLIF vectors.

### Specific Molecules

```bash
easydock_plif -i output.db -o plif_output.txt -d mol_1 mol_2 mol_3
```

## PLIF Similarity

Compare PLIFs to reference interactions:

```bash
easydock_plif -i output.db -o similarity.txt \
    --ref_plif ala31.a.hydrophobic asp86.a.cationic
```

**PLIF naming format:** `[residue_name][residue_number].[chain].[interaction_type]`

**Common interaction types:**

- `hydrophobic`: Hydrophobic interaction
- `hbdonor`: Hydrogen bond donor
- `hbacceptor`: Hydrogen bond acceptor
- `cationic`: Cation-π interaction
- `anionic`: Anion-π interaction


## Combined Operations

Calculate and compare in one command:

```bash
easydock_plif -i output.db -p protein.pdb -o similarity.txt \
    --ref_plif ala31.a.hydrophobic asp86.a.cationic gly45.a.hbdonor -c 4
```

## Workflow Example

**1. Get reference PLIF from known ligand:**
```bash
easydock_plif -i reference.sdf -p protein.pdb -o reference_plif.txt -c 4
```

**2. Identify key interactions from output:**
```
# Example output:
mol_ref: ala31.a.hydrophobic asp86.a.cationic gly45.a.hbdonor
```

**3. Calculate similarity for docked molecules:**
```bash
easydock_plif -i output.db -p protein.pdb -o similarity.txt \
    --ref_plif ala31.a.hydrophobic asp86.a.cationic gly45.a.hbdonor -c 4
```

**4. Extract top scoring by PLIF similarity:**
```bash
get_sdf_from_easydock -i output.db -o output.sdf \
    --fields plif_similarity \
    --add_sql 'plif_similarity > 0.8'
```

## Important Notes

!!! warning "Protein PDB Storage"
    The protein PDB file is stored in the database on first use:
    
    - Used for reference and consistency
    - If you specify a different protein later in the command line, it replaces the stored file
    - **All previously computed PLIFs are erased** without warning
    - This ensures data consistency
    
    **Best practice:** Use a protein file only for the first running of the script `easydock_plif`. ON later runs it will use the protein stored in the database be default.

!!! info "Hydrogen Atoms Required"
    The protein PDB file must include all hydrogen atoms for correct PLIF calculation.

## PLIF Analysis Tips

**Filter by key interactions:**
```bash
# Find molecules with specific interaction
easydock_plif -i output.db -o results.txt -c 4
grep "asp86.a.cationic" results.txt
```

**Combine with docking scores:**
```bash
get_sdf_from_easydock -i output.db -o results.sdf \
    --fields docking_score plif_similarity \
    --add_sql 'docking_score < -8 AND plif_similarity > 0.7'
```
