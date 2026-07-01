#!/usr/bin/env bash
set -euo pipefail

for f in inputs/winners/mol_*/docked.pdbqt; do
  [[ -e "$f" ]] || { echo "no docked.pdbqt found under inputs/winners/mol_*/" >&2; break; }
  d=$(dirname "$f")                       # e.g. inputs/cmpd01
  echo "== $(basename "$d") =="
  ./scripts/pdbqt2sdf.sh "$f" "$d/winner.sdf"
done
