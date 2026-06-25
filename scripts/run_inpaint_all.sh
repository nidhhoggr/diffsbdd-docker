#!/usr/bin/env bash
# run_inpaint_all.sh — run DiffSBDD inpaint over every inputs/cmpd*/ seed.
# Resumable (skips seeds whose output already exists), logs per seed.
#
# Run from the repo root (the bind-mounted workspace). All paths are relative
# to that root because compose mounts .:/workspace.
#
# Env knobs:
#   USE_GPU=1            use the GPU compose override (default: 1)
#   RECEPTOR=path        shared apo receptor PDB (default: inputs/structures/receptor.pdb)
#   N_SAMPLES=50         samples per seed (drop to 25 if you hit CUDA OOM)
#   RESAMPLINGS=20       harmonization steps (quality lever)
#   TIMESTEPS=50         denoising steps
#   ADD_N_NODES=5        new atoms grown onto the fixed hit
#   FORCE=1              re-run seeds even if output exists

set -euo pipefail

USE_GPU="${USE_GPU:-1}"
RECEPTOR="${RECEPTOR:-inputs/structures/receptor.pdb}"
N_SAMPLES="${N_SAMPLES:-50}"
RESAMPLINGS="${RESAMPLINGS:-20}"
TIMESTEPS="${TIMESTEPS:-50}"
ADD_N_NODES="${ADD_N_NODES:-5}"
FORCE="${FORCE:-0}"

[[ -f "$RECEPTOR" ]] || { echo "error: receptor not found: $RECEPTOR" >&2
                          echo "  set RECEPTOR=... to your shared apo receptor PDB" >&2; exit 1; }

# Build the compose invocation (GPU override on/off).
COMPOSE=(docker compose -f docker-compose.yml)
[[ "$USE_GPU" == "1" ]] && COMPOSE+=(-f docker-compose.gpu.yml)

mkdir -p logs
EXTRA="--n_samples $N_SAMPLES --resamplings $RESAMPLINGS --timesteps $TIMESTEPS --sanitize"
echo "receptor=$RECEPTOR  gpu=$USE_GPU  extra='$EXTRA'  add_n_nodes=$ADD_N_NODES"
echo

shopt -s nullglob
ok=0; skip=0; fail=0
for d in inputs/cmpd*/; do
  d="${d%/}"
  name="$(basename "$d")"
  ref="$d/docked.sdf"
  out="$d/elaborated.sdf"

  [[ -s "$ref" ]] || { echo "!! $name: missing $ref — skipping"   >&2; ((fail++)); continue; }
  if [[ -s "$out" && "$FORCE" != "1" ]]; then
    echo "== $name: output exists, skipping (FORCE=1 to redo)"; ((skip++)); continue
  fi

  echo "== $name: generating -> $out"
  if "${COMPOSE[@]}" run --rm \
        -e INPAINT_PDB="$RECEPTOR" \
        -e INPAINT_REF_LIGAND="$ref" \
        -e FIX_ATOMS="$ref" \
        -e INPAINT_OUTFILE="$out" \
        -e CENTER=ligand \
        -e ADD_N_NODES="$ADD_N_NODES" \
        -e INPAINT_EXTRA="$EXTRA" \
        inpaint >"logs/${name}.log" 2>&1; then
    n="$(grep -cF '$$$$' "$out" 2>/dev/null || echo 0)"
    echo "   done: $n molecules"
    ((ok++))
  else
    echo "!! $name: inpaint failed — see logs/${name}.log" >&2
    ((fail++))
  fi
done

echo
echo "summary: $ok generated, $skip skipped, $fail failed"
