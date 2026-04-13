# Client-Server Docking Protocol

This page is intended for **developers** who want to implement a custom docking server compatible with EasyDock's `--program server` mode. See [Usage Guide](usage.md#server-based-docking) for how to configure and run existing servers.

## Overview

EasyDock communicates with a docking server through a persistent subprocess. The server is started once per worker process and remains alive for the entire run. All communication uses **JSON Lines over STDIN/STDOUT**: each message is a JSON object on a single newline-terminated line. The request/response cycle is strictly serialized — only one request is in flight at a time.

Non-JSON stdout lines from the server are silently ignored. Stderr is captured and included in error messages on failure.

---

## Request Message (client → server)

```json
{
  "id": <int>,
  "command": "<string>",
  "payload": { ... }
}
```

| Field | Type | Description |
|---|---|---|
| `id` | int | Auto-incrementing request identifier (starts at 1). Used to match responses. |
| `command` | string | The operation name: `"info"`, `"init"`, or `"dock"`. |
| `payload` | object | Command-specific parameters. Empty object `{}` if none. |

---

## Response Message (server → client)

```json
{
  "id": <int>,
  "payload": {
    "status": "<string>",
    ...
  }
}
```

| Field | Location | Type | Description |
|---|---|---|---|
| `id` | top level | int | Must match the request `id`. Mismatched responses are discarded with a warning. |
| `payload` | top level | object | Result object. Checked first for `status` and `error`. |
| `status` | `payload` | string or bool | Accepted OK values: absent/`null`, `true`, `"ok"`, `"OK"`. Anything else raises `RuntimeError`. |
| `error` | `payload` | string | Error description, included in the exception on failure. |

---

## Command 1 — Info Request

Sent **before** `init` to query the server's default operating parameters. The client uses these defaults to configure ligand preparation and result parsing without requiring the user to specify them manually in the config file.

**Request:**

```json
{"id": 1, "command": "info", "payload": {}}
```

**Response payload fields:**

| Field | Type | Description |
|---|---|---|
| `ligand_in_format` | string | Ligand input format accepted by the server: `"smiles"`, `"mol"`, or `"pdbqt"`. |
| `ligand_out_format` | string | Pose output format returned by the server: `"sdf"` or `"pdbqt"`. |
| `score_mode` | string | Score direction: `"min"` (lower is better, Vina) or `"max"` (higher is better, RTMScore). |
| `batch_size` | int | Preferred number of molecules per dock request. |

**Example response (CarsiDock):**

```json
{
  "id": 1,
  "payload": {
    "ligand_in_format": "smiles",
    "ligand_out_format": "sdf",
    "score_mode": "max",
    "batch_size": 10
  }
}
```

The `info` command must be implemented but may return an empty payload — the client falls back to PDBQT defaults in that case.

---

## Command 2 — Init Request

Sent **once** after the server starts. Command name: `"init"`.

**Payload** is the content of the `init_server:` section from the config YAML — receptor files, pocket definition, and any program-specific parameters.

**Example request (CarsiDock):**

```json
{
  "id": 2,
  "command": "init",
  "payload": {
    "protein": "/path/to/protein.pdb",
    "reflig": "/path/to/reference_ligand.sdf",
    "num_conformer": 5
  }
}
```

### CarsiDock-specific init payload fields

| Field | Type | Default | Description |
|---|---|---|---|
| `protein` | string | required | Protein PDB file path |
| `reflig` | string | required | Reference ligand file (defines binding pocket) |
| `cuda_convert` | bool | `true` | Use CUDA-based LBFGSB optimization via `pydock` |
| `num_threads` | int | `1` | Number of docking threads |
| `num_conformer` | int | `5` | Number of conformers per ligand (minimum 5 enforced internally) |

### Vina-GPU Server-specific init payload fields

| Field | Type | Default | Description |
|---|---|---|---|
| `protein` | string | required | Protein PDBQT file path |
| `protein_setup` | string | required | Grid box definition file |
| `program` | string | `"vina-gpu"` | Program variant: `"vina-gpu"`, `"qvina-gpu"`, or `"qvinaw-gpu"` |
| `n_poses` | int | `9` | Number of output poses |
| `thread` | int | `8000` | Number of OpenCL threads |
| `seed` | int | `0` | Random seed |

**Expected response:**

```json
{"id": 2, "payload": {"status": "ok"}}
```

---

## Command 3 — Dock Request

Sent for each batch of molecules. Command name: `"dock"`.

**Payload** is a flat dict mapping instance name to a single ligand representation string:

```json
{
  "id": 3,
  "command": "dock",
  "payload": {
    "mol_001": "<smiles_or_mol_block_or_pdbqt_string>",
    "mol_002": "<smiles_or_mol_block_or_pdbqt_string>"
  }
}
```

| Key | Type | Description |
|---|---|---|
| Instance name | string | Molecule name from the RDKit mol `_Name` property. For ring-sampled conformers, the name is `<mol_id>__<i>` (e.g. `mol_001__0`). |
| Value | string | Ligand representation. Format determined by `ligand_in_format` from the server's `info` response. |

### Ligand formats

| `ligand_in_format` | Value |
|---|---|
| `"smiles"` | SMILES string |
| `"mol"` | MDL mol block string |
| `"pdbqt"` | PDBQT-format string |

---

## Dock Response

```json
{
  "id": 3,
  "payload": {
    "status": "ok",
    "results": {
      "mol_001": {
        "docking_score": 7.4,
        "mol_block": "mol_001\n...",
        "raw_block": "...\n$$$$\n..."
      },
      "mol_002": null
    }
  }
}
```

The `results` value is a flat dict mapping each instance name back to a result object or `null` on failure:

| Field | Type | Description |
|---|---|---|
| `results` | object | Dict `{instance_name: result_or_null}`.|
| `docking_score` | float | Numeric score. Non-parseable values cause the instance to be skipped. |
| `mol_block` | string | Best-pose MDL mol block. If present and non-empty, used directly; otherwise converted from `raw_block`. |
| `raw_block` | string | Raw pose data — multi-conformer SDF or PDBQT text. Stored in the database for later retrieval. |

Setting an instance's result to `null` (or omitting it) causes the client to log a warning and return `None` for that molecule.

---

## Post-processing on the Client

### Pose conversion (when `mol_block` is absent)

| `ligand_out_format` | Behaviour |
|---|---|
| `"pdbqt"` | `raw_block` must contain `"MODEL"`; converted to mol block via Meeko (`pdbqt2molblock`). |
| `"sdf"` | First SDF record (`raw_block.split('$$$$\n')[0]`) becomes `mol_block`; full `raw_block` kept. |
| other | `NotImplementedError` |

### Score direction

| `score_mode` | Behaviour |
|---|---|
| `"min"` (default) | Lowest score wins — Vina kcal/mol |
| `"max"` | Highest score wins — RTMScore (CarsiDock) |

---

## Error Responses

```json
{"id": 3, "payload": {"status": "error", "error": "No smiles or mol_block provided"}}
```

Unknown commands must return:

```json
{"status": "error", "error": "Unknown command: <name>"}
```

---

## Final Output

`mol_dock` from `server_dock.py` returns a **list** of `(mol_id, output)` tuples — one per input molecule:

```python
(mol_id, {
    "docking_score": 7.4,     # float — best pose score
    "mol_block": "...",        # MDL mol block (first line = mol_id)
    "raw_block": "...",        # raw SDF or PDBQT (if present in response)
    "dock_time": 12.4          # float — average seconds per molecule in the batch
})
```

Returns `(mol_id, None)` for each molecule where ligand preparation fails, the server reports a failure, or no valid poses are returned.

---

## Configuration Reference

The server config YAML has three distinct sections:

### Top-level control keys

Read exclusively by EasyDock; never forwarded to the server.

| Key | Default | Description |
|---|---|---|
| `script_file` | *(required)* | Command to launch the server subprocess. Can be a bare `.sif` path, Docker image name, or a full command string. |
| `server_cwd` | `null` | Working directory for the subprocess. |
| `score_mode` | `"min"` | Score direction: `"min"` or `"max"`. Must be set explicitly when the server uses a non-Vina scoring function (e.g. `max` for RTMScore). |
| `startup_timeout` | `120` | Seconds to wait for the subprocess to become ready. |
| `request_timeout` | `null` | Per-request timeout in seconds (`null` = unlimited). |
| `boron_replacement` | `false` | Replace boron atoms with carbon before docking. |
| `ligand_in_format` | auto | Override the ligand input format detected from the server `info` response. |
| `ligand_out_format` | auto | Override the pose output format detected from the server `info` response. |

### `init_server:` section

Forwarded verbatim as the payload of the `init` command. Contains all server-specific initialization parameters (receptor files, binding pocket definition, program settings, etc.).

```yaml
init_server:
  protein: /path/to/protein.pdb
  reflig: /path/to/reference_ligand.sdf
  num_conformer: 5
```

### `info_server:` section

A dict merged on top of the server's `info` response. Use this to override server defaults without modifying the server itself.

```yaml
info_server:
  batch_size: 5   # reduce batch size from the server's default
```

---

## Implementing a Custom Server

A minimal server must:

1. Read JSON Lines from STDIN in a loop.
2. Respond to `"info"` with its format and scoring defaults.
3. Respond to `"init"` with `{"status": "ok"}` after loading the receptor.
4. Respond to `"dock"` with a `results` dict `{instance_name: pose_dict_or_null}`.
5. Write each response as a JSON Line to STDOUT, wrapping it as `{"id": <id>, "payload": <result>}`.

See `containers/carsidock/carsidock_server.py` and `containers/vinagpu/vinagpu_server.py` in the repository for complete reference implementations.

---

## Container Interface Conventions

EasyDock detects the container type from `script_file` and builds the launch command automatically.

### How EasyDock launches the container

The `script_file` value is parsed as follows:

| `script_file` form | How it is launched |
|---|---|
| Path ending in `.sif` (file exists) | `apptainer run [--nv] [--bind ...] file.sif server` |
| Name that is not an existing file/dir | `docker run -i --rm [--gpus all] [-v ...] image server` |
| Starts with `apptainer` / `singularity` / `docker` | Used as-is (full form — bind mounts and GPU flags are the user's responsibility) |
| Anything else | Run as a plain executable |

`server` is always appended as the first positional argument so the container entrypoint knows to start the JSON Lines server.

### GPU and bind-mount auto-detection

- **GPU**: if `nvidia-smi` is accessible, EasyDock automatically adds `--nv` (Apptainer/Singularity) or `--gpus all` (Docker) to the launch command. No manual configuration is needed.
- **Bind mounts**: EasyDock recursively scans all values in the `init_server:` config section, identifies existing file and directory paths, and mounts their parent directories into the container at the same absolute path. Receptor files and grid definitions passed via `init_server:` are therefore accessible inside the container without any manual `-v` or `--bind` flags.

### Apptainer / Singularity — `%runscript` template

```singularity
%runscript
    case "$1" in
        server)
            shift
            exec python3 /opt/myserver/server.py "$@"
            ;;
        help)
            echo "Usage: apptainer run myserver.sif server"
            ;;
        *)
            exec "$@"
            ;;
    esac
```

### Docker — `ENTRYPOINT` template

```dockerfile
ENTRYPOINT ["/bin/sh", "-c", "\
  case \"$1\" in \
    server) shift; exec python /opt/myserver/server.py \"$@\" ;; \
    help)   echo 'Use: docker run <image> server' ;; \
    *)      exec \"$@\" ;; \
  esac", "--"]
```

!!! note "Stdin must stay open"
    Docker containers must be started with `-i` (`--interactive`) so STDIN remains connected. EasyDock adds `-i` automatically when launching Docker containers in bare form.