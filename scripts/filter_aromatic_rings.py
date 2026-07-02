#!/usr/bin/env python3
"""Filter molecules to no more than N aromatic rings.

Auto-detects input format per file (not per line):
  CSV mode  -- first non-comment line has a comma AND a header containing
               a "smiles" column (e.g. structural_qc.py's
               cmpd,parent_smiles,smiles,core_preserved,tanimoto,scaffold_match
               output). All original columns are preserved and an
               n_arom_rings column is appended, so this can sit directly
               between structural_qc.py and balance_pairs.py without
               dropping parent_smiles/tanimoto/scaffold_match.
  .smi mode -- legacy whitespace-separated, no header:
               "SMILES" | "SMILES  id" | "score  id  SMILES"
Lines that are blank or start with '#' are skipped (legacy mode only).

Usage:
  python filter_aromatic_rings.py in.smi                 # keep <=3, print kept
  python filter_aromatic_rings.py in.smi -n 3 -o kept.smi --rejected rej.smi
  python filter_aromatic_rings.py qc_filtered.csv -n 3 -o dock_me.csv --rejected rej.csv
"""
import sys, csv, argparse
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")

def parse_line(line):
    """Legacy .smi mode: return (smiles, label) from a flexible line format."""
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

def sniff_csv(path):
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            return "," in line and "smiles" in line.lower()
    return False

def run_csv(args):
    rows = list(csv.DictReader(open(args.infile)))
    if not rows or "smiles" not in rows[0]:
        sys.exit(f"error: no 'smiles' column found in {args.infile}")
    fieldnames = list(rows[0].keys()) + ["n_arom_rings"]

    kept, rejected, bad = [], [], 0
    for r in rows:
        m = Chem.MolFromSmiles(r["smiles"])
        if m is None:
            bad += 1
            continue
        n = rdMolDescriptors.CalcNumAromaticRings(m)
        r["n_arom_rings"] = n
        (kept if n <= args.max_aromatic_rings else rejected).append(r)

    out = sys.stdout if args.out == "-" else open(args.out, "w", newline="")
    w = csv.DictWriter(out, fieldnames=fieldnames)
    w.writeheader(); w.writerows(kept)
    if out is not sys.stdout: out.close()
    if args.rejected:
        with open(args.rejected, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader(); w.writerows(rejected)

    sys.stderr.write(
        f"kept {len(kept)} (<= {args.max_aromatic_rings} aromatic rings), "
        f"rejected {len(rejected)}, unparseable {bad}\n")

def run_smi(args):
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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("-n", "--max-aromatic-rings", type=int, default=3)
    ap.add_argument("-o", "--out", default="-", help="kept output (default stdout)")
    ap.add_argument("--rejected", help="optional file for filtered-out molecules")
    ap.add_argument("--format", choices=["auto", "csv", "smi"], default="auto",
                     help="force input format instead of auto-detecting (default: auto)")
    args = ap.parse_args()

    is_csv = (args.format == "csv") or (args.format == "auto" and sniff_csv(args.infile))
    (run_csv if is_csv else run_smi)(args)

if __name__ == "__main__":
    main()
