#!/usr/bin/env python3
"""Assess per-seed balance + pairing quality of the alert-free elaborations,
and preview a balanced TL pair set under a similarity floor and per-seed cap.

Input: structural_qc output (cmpd,parent_smiles,smiles,core_preserved,tanimoto,scaffold_match)

Usage:
  python balance_pairs.py filtered.csv               # report only
  python balance_pairs.py filtered.csv --sim-floor 0.3 --cap 15 \
         --pairs tl_pairs.smi --dock dock_input.smi
"""
import sys, csv, argparse
from collections import defaultdict, Counter
import statistics

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("--sim-floor", type=float, default=0.0,
                    help="min parent->child Tanimoto for a pair to count as a learnable transform")
    ap.add_argument("--cap", type=int, default=0, help="max molecules kept per seed (0 = no cap)")
    ap.add_argument("--require-scaffold", action="store_true",
                    help="keep only scaffold_match=1 (genuine on-series elaborations)")
    ap.add_argument("--pairs"); ap.add_argument("--dock")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.infile)))
    by = defaultdict(list)
    for r in rows:
        r["tanimoto"] = float(r["tanimoto"]); r["scaffold_match"] = int(r["scaffold_match"])
        by[r["cmpd"]].append(r)

    n = len(rows)
    print(f"{'seed':7}{'all':>5}{'scaf':>6}{'medTan':>8}{'share%':>8}")
    scaf_tot = 0
    for s in sorted(by):
        g = by[s]; sm = sum(x["scaffold_match"] for x in g); scaf_tot += sm
        med = statistics.median(x["tanimoto"] for x in g)
        print(f"{s:7}{len(g):5d}{sm:6d}{med:8.2f}{100*len(g)/n:8.1f}")
    top = max(by, key=lambda s: len(by[s]))
    print(f"\nTOTAL {n} mols across {len(by)} seeds | scaffold-retained {scaf_tot} "
          f"({100*scaf_tot/n:.0f}%) | largest seed = {top} ({100*len(by[top])/n:.0f}% of pool, "
          f"{sum(x['scaffold_match'] for x in by[top])} of {scaf_tot} on-scaffold)")

    # preview balanced TL set
    kept = []
    for s, g in by.items():
        pool = [x for x in g if x["tanimoto"] >= a.sim_floor]
        if a.require_scaffold: pool = [x for x in pool if x["scaffold_match"]]
        pool.sort(key=lambda x: -x["tanimoto"])
        if a.cap: pool = pool[:a.cap]
        kept.extend(pool)
    if a.sim_floor or a.cap or a.require_scaffold:
        share = Counter(x["cmpd"] for x in kept)
        topk = max(share, key=share.get) if share else "-"
        print(f"\nafter floor={a.sim_floor} cap={a.cap or '-'} scaffold-only={a.require_scaffold}: "
              f"{len(kept)} pairs, {len(share)} seeds, largest {topk} "
              f"({100*share[topk]/max(len(kept),1):.0f}%)" if share else "  (empty)")

    if a.dock:
        seen=set()
        with open(a.dock,"w") as fh:
            for r in rows:
                if r["smiles"] not in seen:
                    seen.add(r["smiles"]); fh.write(f"{r['smiles']} {r['cmpd']}\n")
        print(f"wrote {len(seen)} unique SMILES -> {a.dock}")
    if a.pairs:
        with open(a.pairs,"w") as fh:
            for r in kept: fh.write(f"{r['parent_smiles']}\t{r['smiles']}\n")
        print(f"wrote {len(kept)} parent->child pairs -> {a.pairs}")

if __name__ == "__main__":
    main()
