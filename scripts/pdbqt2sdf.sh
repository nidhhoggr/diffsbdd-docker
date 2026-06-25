#!/usr/bin/env bash
# pdbqt2sdf.sh — convert a docked AutoDock/Vina pose (.pdbqt) to .sdf inside a
# container, PRESERVING the docked coordinates, for use as a DiffSBDD
# --ref_ligand / --fix_atoms input.
#
# Usage:
#   ./pdbqt2sdf.sh cmpd01.pdbqt                      # -> cmpd01_docked.sdf (top pose)
#   ./pdbqt2sdf.sh cmpd01.pdbqt out/cmpd01_docked.sdf
#
# Env overrides:
#   IMAGE=...    container image that provides `obabel`  (default: htvs-redock:latest)
#   MODEL=all    convert every pose, not just the top-ranked one  (default: 1)
#   KEEP_H=1     keep hydrogens  (default: strip — DiffSBDD full-atom model is heavy-atom)
#
# NOTE: no --gen3d anywhere. obabel reads the explicit pdbqt coordinates as-is,
# which is the whole point: the docked pose IS the pocket location for DiffSBDD.

set -euo pipefail

IN="${1:?usage: $0 <input.pdbqt> [output.sdf]}"
OUT="${2:-${IN%.pdbqt}_docked.sdf}"
IMAGE="${IMAGE:-htvs-redock:latest}"
MODEL="${MODEL:-1}"
KEEP_H="${KEEP_H:-0}"

[[ -f "$IN" ]] || { echo "error: input not found: $IN" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"

# Resolve absolute dirs so we can bind-mount them; pass basenames to the container.
IN_DIR="$(cd "$(dirname "$IN")"  && pwd)"; IN_BASE="$(basename "$IN")"
OUT_DIR="$(cd "$(dirname "$OUT")" && pwd)"; OUT_BASE="$(basename "$OUT")"

# Pose selection: top-ranked model only (-f 1 -l 1) unless MODEL=all.
POSE_FLAGS=(-f 1 -l 1)
[[ "$MODEL" == "all" ]] && POSE_FLAGS=()

# Hydrogen handling: strip by default (-d) for a clean heavy-atom SDF.
H_FLAGS=(-d)
[[ "$KEEP_H" == "1" ]] && H_FLAGS=()

docker run --rm \
  -v "$IN_DIR":/in \
  -v "$OUT_DIR":/out \
  "$IMAGE" \
  obabel "/in/$IN_BASE" -O "/out/$OUT_BASE" "${POSE_FLAGS[@]}" "${H_FLAGS[@]}"

# Sanity check the result so a silent empty file can't slip downstream.
[[ -s "$OUT" ]] || { echo "error: conversion produced no output ($OUT)" >&2; exit 1; }
n_recs="$(grep -cF '$$$$' "$OUT" || true)"
echo "wrote $OUT  (records: ${n_recs:-?})"
[[ "${n_recs:-0}" -gt 1 && "$MODEL" != "all" ]] && \
  echo "warning: >1 record written; expected a single top pose" >&2
exit 0
