# Data Retrieval

## Using --sdf Flag

The simplest way to extract top poses during or after docking:

```bash
easydock -o output.db --sdf
```

Features:
- Returns top scoring pose for each molecule
- Automatically selects best among enumerated stereoisomers
- Creates SDF file of the same name as a database in the same directory
- Includes docking scores in SDF properties

## Using get_sdf_from_easydock

More flexible extraction tool with extensive filtering options.

### Basic Extraction

Extract top poses with docking scores:
```bash
get_sdf_from_easydock -i output.db -o output.sdf --fields docking_score
```

Extract as SMILES:
```bash
get_sdf_from_easydock -i output.db -o output.smi --fields id docking_score
```

### Specific Molecules

By molecule ID:
```bash
get_sdf_from_easydock -i output.db -o output.sdf -d mol_1 mol_4 mol_10
```

From file:
```bash
# Create id_list.txt with one ID per line
get_sdf_from_easydock -i output.db -o output.sdf -d id_list.txt
```

### Specific Poses

Extract second pose:
```bash
get_sdf_from_easydock -i output.db -o output.sdf --poses 2
```

Extract multiple poses:
```bash
get_sdf_from_easydock -i output.db -o output.sdf --poses 1 2 3
```

!!! note
    By default, only top poses are extracted. Specify `--poses` to get others.

### Filtered Extraction

Score threshold:
```bash
get_sdf_from_easydock -i output.db -o output.sdf \
    --fields docking_score \
    --add_sql 'docking_score < -10'
```

### Additional Fields

Extract additional database fields:

```bash
get_sdf_from_easydock -i output.db -o output.sdf --fields docking_score id stereo_id
```

## Combining Options

Extract specific molecules with score filter:

```bash
get_sdf_from_easydock -i output.db -o output.sdf \
    -d mol_1 mol_2 mol_3 \
    --poses 1 2 \
    --fields docking_score \
    --add_sql 'docking_score < -9'
```

## Advanced Flags

### PoseBusters Filter (`--bust`)

Return only poses that passed the PoseBusters check. Bust results must be pre-computed with `easydock_bust` before using this flag (see [Pose Quality Assessment](posebusters.md)).

```bash
get_sdf_from_easydock -i output.db -o output.sdf --bust
```

Combined with explicit poses — each pose is filtered independently:

```bash
get_sdf_from_easydock -i output.db -o output.sdf --poses 1 2 3 --bust
```

### Keep Stereo ID in Names (`-s` / `--keep_stereo_id`)

By default, the stereo ID suffix automatically injected by EasyDock during the preparation step (e.g. `_0`, `_1`, etc) is stripped from molecule names in the output. Use this flag to preserve it:

```bash
get_sdf_from_easydock -i output.db -o output.sdf -s
```

### Debug SQL (`--debug`)

Print the final SQL query to stdout before execution. Useful for troubleshooting custom `--add_sql` expressions:

```bash
get_sdf_from_easydock -i output.db -o output.sdf --add_sql 'docking_score < -10' --debug
```

## Enumerated Stereoisomer Handling

!!! warning "Stereoisomer Limitation"
    `get_sdf_from_easydock` does not automatically select the best pose among enumerated stereoisomers. It treats all database records as independent. 
    
    **Workaround:** If needed, use `easydock -o database.db --sdf` which keep only the best enumerated stereoisomer.
