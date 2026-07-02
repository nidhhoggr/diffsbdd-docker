#!/usr/bin/env python3
"""Filter a SMILES list to molecules with no more than N aromatic rings.

Input formats accepted (auto-detected per line):
  - "SMILES"
  - "SMILES  id"          (whitespace-separated)
  - "score  id  SMILES"   (3-col, like your docking lists)
Lines that are blank or start with '#' are skipped.

Usage:
  python filter_aromatic_rings.py in.smi                 # keep <=3, print kept
  python filter_aromatic_rings.py in.smi -n 3 -o kept.smi --rejected rej.smi
"""
import sys, argparse
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")

def parse_line(line):
    """Return (smiles, label) from a flexible line format."""
    parts = line.split()
    if len(parts) == 1:
        return parts[0], parts[0]
    # 3-col docking style: score id smiles  -> smiles is the token that parses
    if len(parts) >= 3 and Chem.MolFromSmiles(parts[-1]):
        return parts[-1], parts[1]
    # 2-col "smiles id": first token parses as a molecule
    if Chem.MolFromSmiles(parts[0]):
        return parts[0], parts[1] if len(parts) > 1 else parts[0]
    # fallback: last token
    return parts[-1], parts[0]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("-n", "--max-aromatic-rings", type=int, default=3)
    ap.add_argument("-o", "--out", default="-", help="kept SMILES (default stdout)")
    ap.add_argument("--rejected", help="optional file for filtered-out molecules")
    args = ap.parse_args()

    kept, rejected, bad = [], [], 0
    with open(args.infile) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            smi, label = parse_line(line)
            m = Chem.MolFromSmiles(smi)
            if m is None:
                bad += 1
                continue
            n = rdMolDescriptors.CalcNumAromaticRings(m)
            (kept if n <= args.max_aromatic_rings else rejected).append((smi, label, n))

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    for smi, label, n in kept:
        out.write(f"{smi}\t{label}\t{n}\n")
    if out is not sys.stdout:
        out.close()
    if args.rejected:
        with open(args.rejected, "w") as fh:
            for smi, label, n in rejected:
                fh.write(f"{smi}\t{label}\t{n}\n")

    sys.stderr.write(
        f"kept {len(kept)} (<= {args.max_aromatic_rings} aromatic rings), "
        f"rejected {len(rejected)}, unparseable {bad}\n")

if __name__ == "__main__":
    main()
