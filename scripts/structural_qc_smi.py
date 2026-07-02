#!/usr/bin/env python3
"""Structural-alert filter for a SMILES list (no parent column required).

Reads flexible lines: "SMILES", "SMILES id", or "SMILES\tid\tqed".
Drops molecules matching implausible/unstable substructures. Extends the
original alert set with the classes that slipped through before (imide,
hydroxamic acid, azo/azoxy, N-N-O).

Usage: python structural_qc_smi.py in.smi -o clean.smi [--rejected rej.smi]
"""
import sys, argparse
from rdkit import Chem
from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")

ALERTS = {
 "peroxide_OO":     "[OX2]-[OX2]",
 "N-halogen":       "[NX3]-[F,Cl,Br,I]",
 "acyl_halide":     "[CX3](=O)[F,Cl,Br,I]",
 "cumulene/ketene": "[$([CX2](=*)=*)]",
 "allene":          "[#6]=[#6]=[#6]",
 "gem-diol":        "[CX4]([OX2H])[OX2H]",
 "gem-triol_SOOO":  "[SX4]([OX2H])([OX2H])[OX2H]",   # nonsense S(O)(O)O
 "hemiaminal":      "[CX4]([OX2H])[NX3]",
 "hydroxamic/N-OH": "[NX3][OX2H]",                    # now INCLUDED (was excluded before)
 "imide":           "[CX3](=O)[NX3][CX3]=O",
 "azo/azoxy":       "[#7]=[#7]",
 "N-N-O":           "[#7]-[#7]-[OX2]",
}
ALERTS = {k: Chem.MolFromSmarts(v) for k, v in ALERTS.items()}

def parse(line):
    p = line.split()
    if not p: return None, None
    smi = p[0]
    return smi, (p[1] if len(p) > 1 else smi)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile"); ap.add_argument("-o", "--out", default="clean.smi")
    ap.add_argument("--rejected")
    a = ap.parse_args()
    from collections import Counter
    fired = Counter(); kept = []; rej = []; bad = 0
    for line in open(a.infile):
        line = line.strip()
        if not line or line.startswith("#"): continue
        smi, label = parse(line)
        m = Chem.MolFromSmiles(smi)
        if m is None: bad += 1; continue
        hits = [n for n, patt in ALERTS.items() if patt and m.HasSubstructMatch(patt)]
        if hits:
            for h in hits: fired[h] += 1
            rej.append((smi, label, ";".join(hits)))
        else:
            kept.append((smi, label))
    with open(a.out, "w") as fh:
        for smi, label in kept: fh.write(f"{smi}\t{label}\n")
    if a.rejected:
        with open(a.rejected, "w") as fh:
            for smi, label, h in rej: fh.write(f"{smi}\t{label}\t{h}\n")
    sys.stderr.write(f"kept {len(kept)}, rejected {len(rej)}, unparseable {bad}\n")
    sys.stderr.write(f"alerts: {dict(fired.most_common())}\n")

if __name__ == "__main__":
    main()
