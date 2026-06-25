#!/usr/bin/env bash
# Download a DiffSBDD pretrained checkpoint from Zenodo into ./checkpoints.
# Driven by host.env:  CKPT (destination path) and CKPT_URL (source).
# Skips the download if the file already exists.
set -euo pipefail

CKPT="${CKPT:-checkpoints/crossdocked_fullatom_cond.ckpt}"
CKPT_URL="${CKPT_URL:-https://zenodo.org/record/8183747/files/crossdocked_fullatom_cond.ckpt?download=1}"

dest_dir="$(dirname "$CKPT")"
mkdir -p "$dest_dir"

if [[ -s "$CKPT" ]]; then
  echo "[fetch] already present: $CKPT ($(du -h "$CKPT" | cut -f1))"
  exit 0
fi

echo "[fetch] downloading -> $CKPT"
echo "[fetch] from        -> $CKPT_URL"
# -O sets the output name (Zenodo URLs carry a ?download=1 query string).
wget --no-verbose --show-progress -O "$CKPT" "$CKPT_URL"
echo "[fetch] done: $CKPT ($(du -h "$CKPT" | cut -f1))"
