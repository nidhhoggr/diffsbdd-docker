# Inpaint → QC → Docking Prep

Pipeline for running DiffSBDD `inpaint` over a seed batch and getting
the output to a clean, dockable SMILES set. Covers `batch_inpaint.sh`
→ `triage_elaborated.py` → `structural_qc.py` →
`filter_aromatic_rings.py`, with optional `balance_pairs.py` if the
batch is also going to be used as mol2mol TL pairs.

## Directory conventions

Each campaign lives in its own subfolder under `inputs/`, e.g.
`inputs/mol2mol_quinoline2/`:

```
inputs/mol2mol_quinoline2/
├── compounds.smi          # seed manifest (see format below)
├── mol_0000/
│   ├── docked.sdf           # redocked pose for this seed (input to inpaint)
│   └── elaborated.sdf       # raw DiffSBDD inpaint output for this seed
├── mol_0003/
│   ├── docked.sdf
│   └── elaborated.sdf
└── ...
```

- Directory names are `mol_NNNN` (4-digit zero-padded), matching what
  the `redock_vina` script emits.
- `docked.sdf` is the redocked pose for that seed — produced from
  `docked.pdbqt` via `pdbqt2sdf.sh` — and is the direct input to
  `batch_inpaint.sh` (Step 0).
- Each `mol_NNNN/elaborated.sdf` is the raw, unfiltered inpaint output
  for that seed — it will contain invalid records, duplicates, and
  copies identical to the parent; that's expected and is what
  `triage_elaborated.py` cleans up.

### `compounds.smi` format

Whitespace/tab-separated, **no header**, 3 columns:

```
mol_id  docking_score  smiles
```

```
mol_0003        -12.3   O=C1NC(Nc2c1cccc2)c1ccc(cc1)c1cc(cc(c1)C(F)(F)F)C(F)(F)F
mol_0000        -12.1   CNc1nccc(n1)c1cccc(c1)c1cc(cc(c1)C(F)(F)F)C(F)(F)F
```

- `mol_id` must be **byte-identical** to the corresponding directory
  name — matching is by string equality, not position or numeric
  parsing.
- `docking_score` is carried through as metadata only; it isn't parsed
  as a float anywhere in this stage, so `-12` and `-12.0` are both
  fine.
- `smiles` is the **parent** compound that was fed into `inpaint.py`
  for that seed (i.e. the fixed core / re-docked pose that DiffSBDD
  elaborated around).

Blank lines and lines starting with `#` are skipped.

## Step 0 — `batch_inpaint.sh`

**Purpose:** run DiffSBDD's `inpaint.py` (via the `inpaint` Compose
service) over every seed directory, using each seed's own redocked
pose as both the pocket definition and the fixed substructure to
elaborate around.

```
scripts/batch_inpaint.sh inputs/mol2mol_quinoline2
```

Run from the repo root, inside the container (`./scripts/interactive.bash`
drops you into it) or via `docker compose run`.

**Per-seed I/O** (relative to `INPUTS_DIR`, arg 1, default `inputs`):

| | path | role |
|---|---|---|
| in | `<mol_id>/docked.sdf` | `--ref_ligand` (defines the pocket) **and** `--fix_atoms` (the substructure held fixed) |
| out | `<mol_id>/elaborated.sdf` | raw inpaint output — input to Step 1 |

**Shared parameters** come from `host.env` (`INPAINT_PDB` — the
receptor structure, same for every seed in the batch; `CENTER`;
`ADD_N_NODES`; `CKPT`; `INPAINT_EXTRA` for any raw extra flags). Only
the three per-seed paths (`INPAINT_REF_LIGAND`, `FIX_ATOMS`,
`INPAINT_OUTFILE`) are overridden per-directory via `docker compose
run -e`, which takes precedence over `host.env`.

**Resumable:** any `mol_*` whose `elaborated.sdf` already exists
(non-empty) is skipped, so a batch can be safely re-run after a
partial failure. Force a full redo with:

```
FORCE=1 scripts/batch_inpaint.sh inputs/mol2mol_quinoline2
```

Logs go to `logs/inpaint/<mol_id>.log` per seed, plus a running
`logs/inpaint/summary.log` with one OK/SKIP/FAIL line per seed. Exit
status is always 0 — check the printed `Done: N ok, N skipped, N
failed.` line (and `Failed: ...` list) rather than the shell's `$?`.

A seed is marked `FAIL` (not just skipped) if `docked.sdf` is missing
or empty going in, or if the `docker compose run` for that seed exits
non-zero — check that seed's log file for the underlying DiffSBDD
error in the latter case.

> **Note:** the script's own comments and the "no compound dirs"
> error message still say `inputs/cmpd*/` — that's stale wording left
> over from the old convention. The actual glob on disk is `mol*/`,
> which is what matches your directory layout; nothing to change on
> your end.

## Step 1 — `triage_elaborated.py`

**Purpose:** inventory + clean the raw inpaint output. Sanitizes each
record, drops disconnected fragments, dedupes, drops anything
identical to the parent, and flags whether the fixed core survived
elaboration (`core_preserved`).

```
python ./scripts/triage_elaborated.py \
    --inputs inputs/mol2mol_quinoline2/ \
    --compounds inputs/mol2mol_quinoline2/compounds.smi \
    --out inputs/mol2mol_quinoline2/elaborated_clean.csv
```

**Output:** `elaborated_clean.csv` with columns
`cmpd, parent_smiles, smiles, core_preserved`.

Per-seed console summary (`raw / valid / uniq / ≠parent / core✓`) is
worth scanning before moving on — a seed with a `≠parent` count under
~10 is flagged as low-yield in the summary line and may not be worth
carrying forward.

**Note:** lookup is keyed by `mol_id` against `compounds.smi`, not by
line position. A `mol_*` directory with no matching row in
`compounds.smi` is skipped with a `WARN`, not silently mismatched.

## Step 2 — `structural_qc.py`

**Purpose:** chemistry sanity gate on top of the triage pass. Flags
implausible/unstable substructures common in diffusion-model output
(peroxides, N-halogens, acyl halides, cumulenes/allenes, gem-diols,
hemiaminals, N-O amines) that an SA score alone won't reliably catch.
Also scores drift from parent: exact core match, Murcko scaffold
match, and Tanimoto similarity.

```
python ./scripts/structural_qc.py \
    inputs/mol2mol_quinoline2/elaborated_clean.csv \
    -o inputs/mol2mol_quinoline2/qc_filtered.csv
```

**Output:** `qc_filtered.csv` — same columns as the input, plus
`tanimoto` and `scaffold_match`. Only alert-free rows are written.

Console output reports per-seed `core / scaf / clean` counts and which
alert classes fired overall. A seed with `core=0, scaf=0` across all
its rows means every elaboration for that parent drifted to an
entirely different scaffold — not necessarily wrong, but worth a
manual look before docking, since none of those molecules represent
"same series, decorated."

**Variant:** `structural_qc_smi.py` covers a broader alert set (adds
imide, hydroxamic acid, azo/azoxy, N-N-O) but works on a bare SMILES
list with no parent/core-retention scoring. For small batches it can
be worth running as a second pass to catch anything the primary QC
step's narrower alert set missed.

## Step 3 — `filter_aromatic_rings.py`

**Purpose:** cap aromatic ring count before docking, consistent with
the scoring config's `NumAromaticRings` penalty (validated cap: ≤3
rings).

```
python ./scripts/filter_aromatic_rings.py \
    inputs/mol2mol_quinoline2/qc_filtered.csv -n 3 \
    -o inputs/mol2mol_quinoline2/dock_me.csv \
    --rejected inputs/mol2mol_quinoline2/rej_rings.csv
```

**Output:** `dock_me.csv` — same columns as the input CSV, plus
`n_arom_rings`. Rejected rows (over the cap) go to the `--rejected`
file in the same format if specified.

**Note:** this script auto-detects its input format. Point it at a
comma-delimited CSV with a `smiles` header column (i.e. `structural_qc.py`'s
output) and it preserves every column; point it at a legacy
whitespace `.smi` file (`SMILES`, `SMILES id`, or `score id SMILES`,
no header) and it falls back to that format. Force one or the other
with `--format csv` / `--format smi` if needed.

## Optional — `balance_pairs.py`

Only needed if this batch is also going to be used as mol2mol
transfer-learning pairs (parent → elaborated child), not just docked.
Takes `structural_qc.py`'s (or `filter_aromatic_rings.py`'s, since it
preserves the same columns) output and reports per-seed balance —
count, scaffold retention, median Tanimoto, share of the pool — then
previews a balanced pair set under a similarity floor and per-seed
cap.

```
# report only
python ./scripts/balance_pairs.py inputs/mol2mol_quinoline2/dock_me.csv

# emit a capped, floor-filtered pair set + dedupe'd dock list
python ./scripts/balance_pairs.py inputs/mol2mol_quinoline2/dock_me.csv \
    --sim-floor 0.3 --cap 15 \
    --pairs inputs/mol2mol_quinoline2/tl_pairs.smi \
    --dock inputs/mol2mol_quinoline2/dock_input.smi
```

Use `--require-scaffold` to keep only genuine on-series elaborations
(`scaffold_match=1`) if a seed's off-scaffold drift (see Step 2) makes
that a concern for TL pair quality.

## After this

`dock_me.csv` (or `dock_input.smi` from `balance_pairs.py`) is the
SMILES set ready for 3D prep and docking — 2D SMILES → 3D via the
usual Open Babel `--gen3d` pass, then into DockStream/Vina. This is
the reverse direction from `pdbqt2sdf.sh`, which converts docked poses
*back* to SDF post-dock for use as the next round's `--ref_ligand` /
`--fix_atoms` input.
