#!/usr/bin/env bash
# Run DiffSBDD `inpaint` over every compound directory under inputs/cmpd*/.
#
# Per-compound layout expected (relative to the diffsbdd-docker/ project root,
# which is mounted at /workspace):
#   inputs/<cmpd>/docked.sdf   -> used as BOTH --ref_ligand (pocket def) and
#                                 --fix_atoms (fixed substructure to elaborate)
# Output is written next to it:
#   inputs/<cmpd>/elaborated.sdf
#
# Shared params come from host.env (protein via INPAINT_PDB, ADD_N_NODES,
# INPAINT_EXTRA, etc.). Only the three per-compound paths are overridden here
# with `docker compose run -e`, which takes precedence over env_file.
#
# Usage:
#   scripts/batch_inpaint.sh            # iterates inputs/cmpd*/
#   scripts/batch_inpaint.sh inputs2    # custom inputs dir
#
# Resumable: any compound that already has a non-empty elaborated.sdf is
# skipped. To force a redo of everything, run with FORCE=1:
#   FORCE=1 scripts/batch_inpaint.sh
set -uo pipefail

INPUTS_DIR="${1:-inputs}"
LOG_DIR="logs/inpaint"
mkdir -p "$LOG_DIR"

shopt -s nullglob
dirs=( "$INPUTS_DIR"/cmpd*/ )
if (( ${#dirs[@]} == 0 )); then
  echo "No compound dirs under '$INPUTS_DIR' (expected $INPUTS_DIR/cmpd*/)." >&2
  exit 1
fi

ok=0; fail=0; skip=0; failed=()
for d in "${dirs[@]}"; do
  d="${d%/}"                      # strip trailing slash
  cmpd="$(basename "$d")"
  ref_host="$d/docked.sdf"        # host-side path for the existence check
  out_host="$d/elaborated.sdf"    # host-side path for the skip check

  # Already done? Skip so the batch is resumable. Set FORCE=1 to redo.
  if [[ -s "$out_host" && "${FORCE:-0}" != "1" ]]; then
    echo "[$cmpd] SKIP: $out_host already exists" | tee -a "$LOG_DIR/summary.log"
    ((skip++)); continue
  fi

  if [[ ! -s "$ref_host" ]]; then
    echo "[$cmpd] SKIP: missing $ref_host" | tee -a "$LOG_DIR/summary.log"
    ((fail++)); failed+=("$cmpd"); continue
  fi

  # Container-side (/workspace) paths for the overrides.
  ref="/workspace/$INPUTS_DIR/$cmpd/docked.sdf"
  out="/workspace/$INPUTS_DIR/$cmpd/elaborated.sdf"

  echo "[$cmpd] inpaint -> $out"
  if docker compose run --rm \
        -e INPAINT_OUTFILE="$out" \
        -e INPAINT_REF_LIGAND="$ref" \
        -e FIX_ATOMS="$ref" \
        inpaint > "$LOG_DIR/$cmpd.log" 2>&1; then
    echo "[$cmpd] OK" | tee -a "$LOG_DIR/summary.log"
    ((ok++))
  else
    echo "[$cmpd] FAIL (see $LOG_DIR/$cmpd.log)" | tee -a "$LOG_DIR/summary.log"
    ((fail++)); failed+=("$cmpd")
  fi
done

echo "----"
echo "Done: $ok ok, $skip skipped, $fail failed."
(( fail )) && printf 'Failed: %s\n' "${failed[*]}"
exit 0
