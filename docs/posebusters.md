# Pose Quality Assessment

`easydock_bust` runs [PoseBusters](https://github.com/maabuu/posebusters) on docked molecules stored in an EasyDock database and records pass/fail results in a `bust` table. The `fast` protocol is used. Results are stored in the database as a single indicator: True - if all checks were passed, false - if at least one check was failed. These results can then be used to filter output with `get_sdf_from_easydock --bust`.

## Installation

```bash
pip install posebusters
```

## Basic Usage

```bash
easydock_bust -i output.db -p protein.pdb -c 4
```

On the **first run**, the protein structure is stored in the database and reused automatically. You do not need to supply `-p` again on subsequent runs unless you want to replace the protein.

## Arguments

`-i` / `--input` (required)
:   EasyDock SQLite database.

`-p` / `--protein`
:   Protein PDB file with all hydrogen atoms. Required on the first run; stored in the database and reused on subsequent runs. Providing a different file replaces the stored structure and clears all previously computed bust results.

`-d` / `--ids` (default: all docked)
:   Molecule IDs to process. Can be a space-separated list or a path to a text file with one ID per line.

`--poses` (default: `1`)
:   Pose numbers to check, 1-based.

`-c` / `--ncpu` (default: all CPUs)
:   Number of CPUs.

`-o` / `--output` (default: stdout)
:   Output file path.

`--full`
:   Write all individual PoseBusters check columns to the output. When set, all molecules are re-processed from scratch to produce the full report.

## Output Format

Tab-separated file with columns: `id`, `stereo_id`, `pose`, `result`.

```
id      stereo_id  pose  result
mol_1   0          1     True
mol_2   0          1     False
mol_3   1          1     True
```

With `--full`, additional boolean columns for each individual PoseBusters check are appended.

## Examples

### Check all docked molecules (top pose only)

```bash
easydock_bust -i output.db -p protein.pdb -c 4
```

### Check multiple poses

```bash
easydock_bust -i output.db -p protein.pdb --poses 1 2 3 -c 4
```

### Write full PoseBusters report to file

```bash
easydock_bust -i output.db --poses 1 -c 4 --full -o bust_report.tsv
```

### Process specific molecules

```bash
easydock_bust -i output.db -p protein.pdb -d mol_1 mol_2 mol_3 -c 4
```

Or from a file:

```bash
easydock_bust -i output.db -p protein.pdb -d id_list.txt -c 4
```

### Resumable runs

Already-processed molecules are skipped automatically (unless `--full` is used). You can re-run the same command to add new molecules without reprocessing existing ones:

```bash
# First run: processes all docked molecules
easydock_bust -i output.db -p protein.pdb -c 4

# Later: only new molecules are processed
easydock_bust -i output.db -c 4
```

## Filtering Output with `get_sdf_from_easydock`

After running `easydock_bust`, use the `--bust` flag to keep only poses that passed:

```bash
get_sdf_from_easydock -i output.db -o output.sdf --bust
```

For specific poses:

```bash
get_sdf_from_easydock -i output.db -o output.sdf --poses 1 2 --bust
```

See [Data Retrieval](data_retrieval.md) for more filtering options.

!!! warning "Protein PDB Storage"
    The protein PDB file is stored in the database on first use. If you supply a **different** protein on a later run, all previously computed bust results are erased and the new protein is stored.

    **Best practice:** Supply `-p` only on the first run. On subsequent runs omit it; the stored protein is reused automatically.

!!! info "Hydrogen Atoms Required"
    The protein PDB file must include all hydrogen atoms for correct PoseBusters validation.
