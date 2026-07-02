#!/usr/bin/env python3
"""Annotate/cluster the docking-hit seeds for provenance tracking.

Input : inputs/compounds.smi -- whitespace/tab-separated, no header:
            mol_id  docking_score  smiles
        e.g. "mol_0001  -11.2  CC(=O)Nc1ccc(cc1)..."
        mol_id is used directly as the 'cmpd' key (matches the mol_*
        directory names produced by the redock_vina script), so this file
        can be freely reordered or subset.
Output: a CSV (default stdout) with one row per compound:
        cmpd, score, smiles, murcko, generic_scaffold, mw, clogp,
        n_rings, n_arom_rings, has_CF3, n_halogen, has_amide, biaryl, cluster

Usage:
  python classify_seeds.py [compounds.smi] [-o seeds.csv] [--sim 0.4]
"""
import sys, csv, argparse
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, DataStructs
from rdkit.ML.Cluster import Butina

# Morgan fingerprint: prefer the new generator API, fall back to legacy
# (older RDKit builds lack rdFingerprintGenerator.GetMorganGenerator).
try:
    from rdkit.Chem import rdFingerprintGenerator
    _MFP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    def morgan(m):
        return _MFP.GetFingerprint(m)
except (ImportError, AttributeError):
    def morgan(m):
        return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)

CF3   = Chem.MolFromSmarts("[CX4](F)(F)F")
AMIDE = Chem.MolFromSmarts("[NX3][CX3](=O)")
HALO  = Chem.MolFromSmarts("[F,Cl,Br,I]")

def biaryl(m):  # two aromatic ring atoms directly single-bonded across rings
    ri = m.GetRingInfo()
    aromatic_atoms = {a.GetIdx() for a in m.GetAtoms() if a.GetIsAromatic()}
    for b in m.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if (a1.GetIdx() in aromatic_atoms and a2.GetIdx() in aromatic_atoms
                and not b.IsInRing() and b.GetBondType() == Chem.BondType.SINGLE):
            return True
    return False

def parse(path):
    """mol_id, score, smiles) triples. Requires all 3 columns."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(
                    f"expected 'mol_id  docking_score  smiles', got: {line!r}")
            mol_id, score, smi = parts[0], parts[1], parts[2]
            rows.append((mol_id, score, smi))
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile", nargs="?", default="inputs/compounds.smi")
    ap.add_argument("-o", "--out", default="-")
    ap.add_argument("--sim", type=float, default=0.4,
                    help="Tanimoto similarity cutoff for Butina clustering")
    args = ap.parse_args()

    rows = parse(args.infile)
    mols, recs = [], []
    fps = []
    for mol_id, score, smi in rows:
        m = Chem.MolFromSmiles(smi)
        cmpd = mol_id
        if m is None:
            sys.stderr.write(f"WARN {cmpd}: unparseable SMILES: {smi}\n")
            recs.append(dict(cmpd=cmpd, score=score, smiles=smi)); fps.append(None)
            continue
        murcko  = MurckoScaffold.GetScaffoldForMol(m)
        generic = MurckoScaffold.MakeScaffoldGeneric(murcko)
        recs.append(dict(
            cmpd=cmpd, score=score, smiles=Chem.MolToSmiles(m),
            murcko=Chem.MolToSmiles(murcko),
            generic_scaffold=Chem.MolToSmiles(generic),
            mw=round(Descriptors.MolWt(m), 1),
            clogp=round(Crippen.MolLogP(m), 2),
            n_rings=rdMolDescriptors.CalcNumRings(m),
            n_arom_rings=rdMolDescriptors.CalcNumAromaticRings(m),
            has_CF3=int(m.HasSubstructMatch(CF3)),
            n_halogen=len(m.GetSubstructMatches(HALO)),
            has_amide=int(m.HasSubstructMatch(AMIDE)),
            biaryl=int(biaryl(m)),
        ))
        fps.append(morgan(m))

    # Butina cluster on valid fps; assign cluster ids back in original order
    valid_idx = [i for i, f in enumerate(fps) if f is not None]
    vfps = [fps[i] for i in valid_idx]
    n = len(vfps)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(vfps[i], vfps[:i])
        dists.extend(1 - s for s in sims)
    clusters = Butina.ClusterData(dists, n, 1 - args.sim, isDistData=True) if n else []
    cl_of = {}
    for cid, members in enumerate(clusters, 1):
        for local in members:
            cl_of[valid_idx[local]] = cid
    for i, r in enumerate(recs):
        r["cluster"] = cl_of.get(i, "")

    cols = ["cmpd","score","cluster","murcko","generic_scaffold","mw","clogp",
            "n_rings","n_arom_rings","has_CF3","n_halogen","has_amide","biaryl","smiles"]
    out = sys.stdout if args.out == "-" else open(args.out, "w", newline="")
    w = csv.DictWriter(out, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    for r in recs: w.writerow(r)
    if out is not sys.stdout: out.close()

    sizes = sorted((len(c) for c in clusters), reverse=True)
    sys.stderr.write(
        f"\n{len(recs)} compounds | {len({r.get('murcko') for r in recs})} distinct Murcko "
        f"scaffolds | {len(clusters)} clusters @ sim {args.sim} (sizes {sizes})\n")

if __name__ == "__main__":
    main()
