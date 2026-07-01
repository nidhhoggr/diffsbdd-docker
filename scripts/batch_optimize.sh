#!/usr/bin/env bash
# Run DiffSBDD `optimize` (evolutionary QED/SA) over selected molecule dirs.
#
# Per-dir layout expected (relative to the diffsbdd-docker/ project root,
# mounted at /workspace):
#   <INPUTS_DIR>/<mol>/winner.sdf   -> the docked pose to optimize (--ref_ligand)
# Output is written next to it:
#   <INPUTS_DIR>/<mol>/optimized.sdf
#
# Shared params come from host.env (CKPT, OBJECTIVE, POPULATION_SIZE,
# EVOLUTION_STEPS, TOP_K, TIMESTEPS, OPT_EXTRA). Only the per-molecule
# receptor/ref/out paths are overridden here with `docker compose run -e`.
#
# Usage:
#   scripts/batch_optimize.sh                  # iterates winners/mol*/
#   scripts/batch_optimize.sh inputs/winners   # custom dir
#   OBJECTIVE=sa scripts/batch_optimize.sh     # override objective for the run
#   FORCE=1 scripts/batch_optimize.sh          # redo even if optimized.sdf exists
#
# Receptor is shared across all molecules (same pocket you docked into):
RECEPTOR="${RECEPTOR:-inputs/structures/7KEW_SM_GP_protomer.pdb}"

set -uo pipefail
INPUTS_DIR="${1:-inputs/winners}"
LOG_DIR="logs/optimize"
mkdir -p "$LOG_DIR"

[[ -f "$RECEPTOR" ]] || { echo "receptor not found: $RECEPTOR (set RECEPTOR=...)" >&2; exit 1; }

shopt -s nullglob
dirs=( "$INPUTS_DIR"/mol*/ )
if (( ${#dirs[@]} == 0 )); then
  echo "No molecule dirs under '$INPUTS_DIR' (expected $INPUTS_DIR/mol*/)." >&2
  exit 1
fi

ok=0; fail=0; skip=0; failed=()
for d in "${dirs[@]}"; do
  d="${d%/}"
  mol="$(basename "$d")"
  ref_host="$d/winner.sdf"
  out_host="$d/optimized.sdf"

  if [[ -s "$out_host" && "${FORCE:-0}" != "1" ]]; then
    echo "[$mol] SKIP: $out_host exists" | tee -a "$LOG_DIR/summary.log"; ((skip++)); continue
  fi
  if [[ ! -s "$ref_host" ]]; then
    echo "[$mol] SKIP: missing $ref_host" | tee -a "$LOG_DIR/summary.log"; ((fail++)); failed+=("$mol"); continue
  fi

  ref="/workspace/$INPUTS_DIR/$mol/winner.sdf"
  out="/workspace/$INPUTS_DIR/$mol/optimized.sdf"
  pdb="/workspace/$RECEPTOR"
  echo "[$mol] optimize ($OBJECTIVE) -> $out"
  if docker compose run --rm \
        -e OPT_PDB="$pdb" \
        -e OPT_REF_LIGAND="$ref" \
        -e OPT_OUTFILE="$out" \
        optimize > "$LOG_DIR/$mol.log" 2>&1; then
    echo "[$mol] OK" | tee -a "$LOG_DIR/summary.log"; ((ok++))
  else
    echo "[$mol] FAIL (see $LOG_DIR/$mol.log)" | tee -a "$LOG_DIR/summary.log"; ((fail++)); failed+=("$mol")
  fi
done

echo "----"
echo "Done: $ok ok, $skip skipped, $fail failed."
(( fail )) && printf 'Failed: %s\n' "${failed[*]}"
exit 0
