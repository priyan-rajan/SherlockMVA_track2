#!/usr/bin/env python3
"""Parse fpocket output for a PDB and emit a ranked TSV of pockets.

fpocket writes:
  <stem>_out/<stem>_info.txt          <- all pockets: Score, Druggability Score, Volume
  <stem>_out/pockets/pocketN_atm.pdb  <- lining atoms per pocket
We rank by Druggability Score and flag whether a residue of interest (ROI, default
1002 = BUB1B p.Asn1002Lys) lines any pocket, so WT vs mutant can be compared.
"""
import sys, os, re, glob

pdb = sys.argv[1]
label = sys.argv[2]
out_tsv = sys.argv[3]
roi = int(sys.argv[4]) if len(sys.argv) > 4 else 1002

base = os.path.basename(pdb)
stem = base[:-4] if base.lower().endswith(".pdb") else base
outdir = os.path.join(os.path.dirname(pdb), stem + "_out")
info = os.path.join(outdir, f"{stem}_info.txt")
pockets_dir = os.path.join(outdir, "pockets")

if not os.path.exists(info):
    sys.stderr.write(f"[fpocket_summary] MISSING {info}\n")
    open(out_tsv, "w").write("protein\tpocket\tdruggability\tvolume\tfpocket_score\tcontains_ROI\n")
    sys.exit(1)

text = open(info).read()
# split into per-pocket blocks
blocks = re.split(r"(?m)^Pocket\s+\d+\s*:\s*$", text)[1:]

def find(pat, s):
    m = re.search(pat, s, re.M)
    return float(m.group(1)) if m else float("nan")

def substr_int(ln):
    try:
        return int(ln[22:26])
    except ValueError:
        return -1

rows = []
for i, b in enumerate(blocks, start=1):
    drugg = find(r"(?m)^\s*Druggability Score\s*:\s*([-\d.]+)", b)
    volume = find(r"(?m)^\s*Volume\s*:\s*([\d.]+)", b)
    score = find(r"(?m)^\s*Score\s*:\s*([-\d.]+)", b)
    atm = os.path.join(pockets_dir, f"pocket{i}_atm.pdb")
    has_roi = "NA"
    if os.path.exists(atm):
        hit = any(substr_int(ln) == roi for ln in open(atm) if ln.startswith(("ATOM", "HETATM")))
        has_roi = "yes" if hit else "no"
    rows.append((label, i, drugg, volume, score, has_roi))

rows.sort(key=lambda r: (r[2] if r[2] == r[2] else -1), reverse=True)
with open(out_tsv, "w") as out:
    out.write("protein\tpocket\tdruggability\tvolume\tfpocket_score\tcontains_ROI\n")
    for r in rows:
        out.write(f"{r[0]}\t{r[1]}\t{r[2]:.3f}\t{r[3]:.1f}\t{r[4]:.3f}\t{r[5]}\n")
n_roi = sum(1 for r in rows if r[5] == "yes")
print(f"[fpocket_summary] {label}: {len(rows)} pockets, {n_roi} line residue {roi} -> {out_tsv}")
