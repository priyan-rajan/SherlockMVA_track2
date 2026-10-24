#!/usr/bin/env python3
"""Fetch an AlphaFold DB model for a UniProt accession (with User-Agent; the
static file endpoint returns 403 without one). Writes <acc>.pdb and a per-residue
pLDDT table <acc>_plddt.tsv."""
import sys, json, requests

acc = sys.argv[1]
out_pdb = sys.argv[2]
out_plddt = sys.argv[3]
UA = {"User-Agent": "Mozilla/5.0 (bub1b-virtual-docking)"}

api = f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"
r = requests.get(api, headers=UA, timeout=60)
r.raise_for_status()
entries = r.json()
entry = entries[0]
pdb_url = entry["pdbUrl"]
print(f"[fetch_af] {acc} entry={entry['entryId']} latest={entry.get('latestVersion')}")
print(f"[fetch_af] downloading {pdb_url}")
pdb = requests.get(pdb_url, headers=UA, timeout=120)
pdb.raise_for_status()
with open(out_pdb, "wb") as fh:
    fh.write(pdb.content)

# per-residue pLDDT from CA B-factors
with open(out_plddt, "w") as out:
    out.write("residue\tchain\tplddt\n")
    for ln in pdb.text.splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
            resseq = int(ln[22:26])
            chain = ln[21]
            b = float(ln[60:66])
            out.write(f"{resseq}\t{chain}\t{b:.2f}\n")
print(f"[fetch_af] wrote {out_pdb} and {out_plddt}")
