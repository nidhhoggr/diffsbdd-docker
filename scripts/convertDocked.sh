#!/usr/bin/env bash
set -euo pipefail

DIR=${1:-}

for f in $DIR/mol*/docked.pdbqt; do
  [[ -e "$f" ]] || { echo "no docked.pdbqt found under $DIR/mol*/" >&2; break; }
  d=$(dirname "$f")                       # e.g. inputs/cmpd01
  echo "== $(basename "$d") =="
  ./scripts/pdbqt2sdf.sh "$f" "$d/docked.sdf"
done
