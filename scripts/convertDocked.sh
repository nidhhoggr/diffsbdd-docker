#!/usr/bin/env bash
set -euo pipefail

for f in inputs/cmpd*/docked.pdbqt; do
  [[ -e "$f" ]] || { echo "no docked.pdbqt found under inputs/cmpd*/" >&2; break; }
  d=$(dirname "$f")                       # e.g. inputs/cmpd01
  echo "== $(basename "$d") =="
  ./scripts/pdbqt2sdf.sh "$f" "$d/docked.sdf"
done
