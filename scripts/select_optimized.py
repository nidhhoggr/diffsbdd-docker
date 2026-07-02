#!/usr/bin/env python3
"""Collect dockable analogs from DiffSBDD optimize CSVs across winners.

For each inputs/.../<mol>/optimized.csv it keeps rows that are:
  - fate == 'survived'          (drops the 'initial' seed row)
  - QED strictly greater than that winner's initial QED (--min-gain to raise the bar)
  - chemically sane (RDKit-parseable, sanitizable, single fragment)
Emits a SMILES file (canonical) + provenance, ready for structural_qc/aromatic
filtering and then docking. De-dupes canonical SMILES across all winners.

Usage:
  python select_optimized.py 'inputs/enamine/winners/*/optimized.csv' -o dock_me.smi
  python select_optimized.py '.../*/optimized.csv' --min-gain 0.05 --top 5 -o dock_me.smi
"""
import sys, csv, glob, argparse
from rdkit import Chem
from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")

def sane(smi):
    m = Chem.MolFromSmiles(smi)          # MolFromSmiles already sanitizes
    if m is None: return None
    c = Chem.MolToSmiles(m)
    return None if "." in c else c        # drop disconnected fragments

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("glob", help="glob for optimized.csv files (quote it)")
    ap.add_argument("-o", "--out", default="dock_me.smi")
    ap.add_argument("--min-gain", type=float, default=0.0,
                    help="min QED improvement over each winner's initial (default: >0)")
    ap.add_argument("--top", type=int, default=0,
                    help="keep only the top-N by QED per winner (0 = all that qualify)")
    a = ap.parse_args()

    files = sorted(glob.glob(a.glob))
    if not files:
        sys.exit(f"no files matched: {a.glob}")
    seen, kept, stats = set(), [], []
    for f in files:
        mol_id = f.split("/")[-2]
        rows = list(csv.DictReader(open(f)))
        init = next((r for r in rows if r.get("fate") == "initial"), None)
        base = float(init["score"]) if init else 0.0
        cand = []
        for r in rows:
            if r.get("fate") != "survived": continue
            q = float(r["score"])
            if q <= base + a.min_gain: continue
            c = sane(r["smiles"])
            if not c: continue
            cand.append((q, c))
        cand.sort(reverse=True)
        if a.top: cand = cand[:a.top]
        n_new = 0
        for q, c in cand:
            if c in seen: continue
            seen.add(c); kept.append((c, mol_id, round(q, 3))); n_new += 1
        stats.append((mol_id, round(base, 3), len(cand), n_new))

    with open(a.out, "w") as fh:
        for c, mol_id, q in kept:
            fh.write(f"{c}\t{mol_id}\t{q}\n")

    print(f"{'winner':12}{'initQED':>8}{'qualify':>8}{'new':>5}")
    for s in stats: print(f"{s[0]:12}{s[1]:8.2f}{s[2]:8d}{s[3]:5d}")
    print(f"\nwrote {len(kept)} unique analogs -> {a.out}  "
          f"(min gain over initial: {a.min_gain}, top/winner: {a.top or 'all'})")

if __name__ == "__main__":
    main()
