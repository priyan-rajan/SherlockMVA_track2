#!/usr/bin/env python3
"""Full survey of every fpocket pocket on a structure, annotated with the mean
pLDDT of its lining residues and the UniProt domain it falls in. Intended to
answer: is there any well-folded, sufficiently druggable cavity that a small
molecule could occupy to *stabilise* the protein?

Usage: pocket_survey.py <pdb> <plddt_tsv> <outdir>
Writes <outdir>/pockets_annotated.tsv (one row per pocket, ranked by druggability).
"""
import sys, os, re, glob, subprocess, collections

pdb, plddt_tsv, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(outdir, exist_ok=True)

# --- pLDDT map (residue -> CA pLDDT) ---
plddt = {}
with open(plddt_tsv) as fh:
    next(fh)
    for ln in fh:
        r, ch, v = ln.split("\t")
        plddt[int(r)] = float(v)

# --- run fpocket ---
subprocess.run(["fpocket", "-f", os.path.basename(pdb)], check=True,
               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
stem = os.path.basename(pdb)[:-4] if pdb.lower().endswith(".pdb") else os.path.basename(pdb)
base = os.path.join(os.path.dirname(pdb), stem + "_out")
info = os.path.join(base, f"{stem}_info.txt")
pockets_dir = os.path.join(base, "pockets")

text = open(info).read()
blocks = re.split(r"(?m)^Pocket\s+\d+\s*:\s*$", text)[1:]

def find(pat, s):
    m = re.search(pat, s, re.M)
    return float(m.group(1)) if m else float("nan")

# --- domain map (UniProt O60566) ---
DOMAINS = [(62, 226, "BUB1_Nterm"), (766, 1050, "Protein_kinase")]
REGIONS = [(111, 118, "NLS"), (152, 185, "KNL1_bind"), (224, 232, "D-box"),
           (579, 580, "caspase_site"), (610, 611, "caspase_site"),
           (772, 780, "ATP_bind"), (795, 795, "ATP_bind"), (882, 882, "active_site")]

def annotate(s, e):
    labs = [n for a, b, n in DOMAINS if s <= b and e >= a]
    labd = [n for a, b, n in REGIONS if s <= b and e >= a]
    return ";".join(labs) if labs else "other", ";".join(labd) if labd else ""

rows = []
for i, b in enumerate(blocks, start=1):
    drugg = find(r"(?m)^\s*Druggability Score\s*:\s*([-\d.]+)", b)
    volume = find(r"(?m)^\s*Volume\s*:\s*([\d.]+)", b)
    score = find(r"(?m)^\s*Score\s*:\s*([-\d.]+)", b)
    atm = os.path.join(pockets_dir, f"pocket{i}_atm.pdb")
    res = []
    if os.path.exists(atm):
        for ln in open(atm):
            if ln.startswith(("ATOM", "HETATM")):
                try:
                    res.append(int(ln[22:26]))
                except ValueError:
                    pass
    res = sorted(set(res))
    if not res:
        continue
    rmin, rmax = res[0], res[-1]
    vals = [plddt[r] for r in res if r in plddt]
    pmean = sum(vals) / len(vals) if vals else float("nan")
    pmin = min(vals) if vals else float("nan")
    fge70 = sum(1 for v in vals if v >= 70) / len(vals) if vals else float("nan")
    dom, reg = annotate(rmin, rmax)
    rows.append(dict(pocket=i, druggability=drugg, volume=volume, score=score,
                     n_res=len(res), res_min=rmin, res_max=rmax,
                     plddt_mean=pmean, plddt_min=pmin, frac_plddt_ge70=fge70,
                     domain=dom, features=reg))

rows.sort(key=lambda r: (r["druggability"] if r["druggability"] == r["druggability"] else -1),
          reverse=True)
cols = ["pocket", "druggability", "volume", "score", "n_res", "res_min", "res_max",
        "plddt_mean", "plddt_min", "frac_plddt_ge70", "domain", "features"]
out = os.path.join(outdir, "pockets_annotated.tsv")
with open(out, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(
            f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
print(f"[pocket_survey] {len(rows)} pockets -> {out}")
