#!/usr/bin/env python3
"""Inventory + clean DiffSBDD inpaint output across inputs/mol_*/elaborated.sdf.

For each seed it reports raw records, how many survive sanitization, how many
are unique, how many actually differ from the parent, and how many still
contain the parent as a substructure (i.e. the fixed core was preserved).
It writes a cleaned, deduped, provenance-tagged SMILES set ready for the SA
gate + re-docking.

inputs/compounds.smi format (whitespace/tab-separated, no header):
    mol_id  docking_score  smiles
e.g.
    mol_0001  -11.2  CC(=O)Nc1ccc(cc1)C(=O)N2CCOCC2

Lookup is keyed by mol_id (matched against the directory name under
--inputs), NOT by line position -- reordering or subsetting compounds.smi
is safe.

Usage:
  python triage_elaborated.py [--inputs inputs] [--compounds inputs/compounds.smi]
                              [--out inputs/elaborated_clean.csv] [--prefix mol_]
"""
import os, csv, glob, argparse
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")          # silence per-molecule sanitize spew

def load_compounds(path):
    """mol_id -> (canonical_smiles, mol). Ignores the docking-score column."""
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(
                    f"expected 'mol_id  docking_score  smiles', got: {line!r}")
            mol_id, _score, smi = parts[0], parts[1], parts[2]
            m = Chem.MolFromSmiles(smi)
            out[mol_id] = (Chem.MolToSmiles(m) if m else None, m)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default="inputs")
    ap.add_argument("--compounds", default="inputs/compounds.smi")
    ap.add_argument("--out", default="inputs/elaborated_clean.csv")
    ap.add_argument("--prefix", default="mol_")
    args = ap.parse_args()

    compounds = load_compounds(args.compounds)

    rows = []          # output: cmpd, parent_smiles, smiles
    print(f"{'seed':10} {'raw':>4} {'valid':>5} {'uniq':>4} {'≠parent':>7} {'core✓':>5}")
    tot = dict(raw=0, valid=0, uniq=0, diff=0, core=0)
    low = []
    for d in sorted(glob.glob(os.path.join(args.inputs, f"{args.prefix}*"))):
        sdf = os.path.join(d, "elaborated.sdf")
        if not os.path.isfile(sdf): continue
        name = os.path.basename(d)
        if name not in compounds:
            print(f"WARN {name}: no matching mol_id in {args.compounds}, skipping")
            continue
        pcanon, pmol = compounds[name]

        raw = valid = 0
        seen = set(); kept = []; core_ok = 0
        for m in Chem.SDMolSupplier(sdf, sanitize=False, removeHs=True):
            raw += 1
            if m is None: continue
            try: Chem.SanitizeMol(m)
            except Exception: continue
            valid += 1
            smi = Chem.MolToSmiles(m)
            if "." in smi: continue                  # disconnected -> drop
            if smi in seen: continue                 # dedupe
            seen.add(smi)
            if pcanon and smi == pcanon: continue    # identical to parent -> drop
            preserved = bool(pmol) and m.HasSubstructMatch(pmol)
            if preserved: core_ok += 1
            kept.append((smi, preserved))

        uniq = len(seen)
        diff = len(kept)
        print(f"{name:10} {raw:4d} {valid:5d} {uniq:4d} {diff:7d} {core_ok:5d}")
        for smi, preserved in kept:
            rows.append(dict(cmpd=name, parent_smiles=pcanon or "", smiles=smi,
                             core_preserved=int(preserved)))
        tot["raw"]+=raw; tot["valid"]+=valid; tot["uniq"]+=uniq
        tot["diff"]+=diff; tot["core"]+=core_ok
        if diff < 10: low.append((name, diff))

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["cmpd","parent_smiles","smiles","core_preserved"])
        w.writeheader(); w.writerows(rows)

    print(f"\nTOTAL   raw={tot['raw']} valid={tot['valid']} unique={tot['uniq']} "
          f"distinct-from-parent={tot['diff']} core-preserved={tot['core']}")
    print(f"wrote {len(rows)} cleaned molecules -> {args.out}")
    if low:
        print("low-yield seeds (<10 usable):", ", ".join(f"{n}({c})" for n,c in low))

if __name__ == "__main__":
    main()
