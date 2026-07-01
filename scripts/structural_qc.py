#!/usr/bin/env python3
"""Structural sanity + core-retention QC for DiffSBDD elaborated output.

Beyond the triage pass (valid/unique/≠parent), this flags chemically
implausible groups that diffusion models emit and that an SA score alone
won't reliably catch, and it measures how far each output drifted from its
parent (exact-core match, scaffold match, Tanimoto).

Usage: python structural_qc.py elaborated_clean.csv [-o filtered.csv]
"""
import sys, csv, argparse
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

# Implausible / unstable substructures common in 3D-diffusion output.
ALERTS = {
 "peroxide_OO":      "[OX2]-[OX2]",
 "N-halogen":        "[NX3]-[F,Cl,Br,I]",
 "acyl_halide":      "[CX3](=O)[F,Cl,Br,I]",
 "cumulene/ketene":  "[$([CX2](=*)=*)]",
 "allene":           "[#6]=[#6]=[#6]",
 "gem-diol":         "[CX4]([OX2H])[OX2H]",
 "hemiaminal":       "[CX4]([OX2H])[NX3]",
 "N-O_amine":        "[NX3;!$([N+](=O)[O-]);!$(N-C=O);!$(N=*)]-[OX2;H1,H0;!$(O=*)]",
}
ALERTS = {k: Chem.MolFromSmarts(v) for k, v in ALERTS.items()}

def morgan(m):
    try:
        from rdkit.Chem import rdFingerprintGenerator as G
        return G.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(m)
    except Exception:
        return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("-o", "--out", default="elaborated_qc.csv")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.infile)))
    per = defaultdict(lambda: dict(n=0, core=0, scaf=0, clean=0, sims=[]))
    alert_tot = defaultdict(int)
    out = []
    pcache = {}

    for r in rows:
        psmi, smi = r["parent_smiles"], r["smiles"]
        seed = r["cmpd"]
        m = Chem.MolFromSmiles(smi)
        if m is None: continue
        if psmi not in pcache:
            pm = Chem.MolFromSmiles(psmi)
            pscaf = MurckoScaffold.GetScaffoldForMol(pm) if pm else None
            pcache[psmi] = (pm, pscaf, morgan(pm) if pm else None)
        pm, pscaf, pfp = pcache[psmi]

        d = per[seed]; d["n"] += 1
        core_exact = bool(pm) and m.HasSubstructMatch(pm)
        scaf_ok    = bool(pscaf) and m.HasSubstructMatch(pscaf)
        sim = DataStructs.TanimotoSimilarity(pfp, morgan(m)) if pfp else 0.0
        d["core"] += core_exact; d["scaf"] += scaf_ok; d["sims"].append(sim)

        fired = [name for name, patt in ALERTS.items() if patt and m.HasSubstructMatch(patt)]
        for f in fired: alert_tot[f] += 1
        if not fired:
            d["clean"] += 1
            out.append({**r, "tanimoto": round(sim,3), "scaffold_match": int(scaf_ok)})

    import statistics
    print(f"{'seed':7}{'n':>4}{'core':>5}{'scaf':>5}{'clean':>6}{'medTan':>8}")
    tot = dict(n=0, core=0, scaf=0, clean=0)
    for seed in sorted(per):
        d = per[seed]; med = statistics.median(d["sims"]) if d["sims"] else 0
        print(f"{seed:7}{d['n']:4d}{d['core']:5d}{d['scaf']:5d}{d['clean']:6d}{med:8.2f}")
        for k in tot: tot[k]+=d[k]
    print(f"\nTOTAL n={tot['n']}  core-exact={tot['core']}  scaffold-match={tot['scaf']}  "
          f"alert-free={tot['clean']}  ({100*tot['clean']/max(tot['n'],1):.0f}% clean)")
    print("alerts fired:", dict(sorted(alert_tot.items(), key=lambda x:-x[1])))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()) if out else
                           ["cmpd","parent_smiles","smiles","core_preserved","tanimoto","scaffold_match"])
        w.writeheader(); w.writerows(out)
    print(f"\nwrote {len(out)} alert-free molecules -> {args.out}")

if __name__ == "__main__":
    main()
