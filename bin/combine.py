#!/usr/bin/env python3
"""Concatenate per-region fpocket summary TSVs into one ranked table."""
import sys

out_path = sys.argv[1]
inputs = sys.argv[2:]
header = None
rows = []
for f in inputs:
    with open(f) as fh:
        lines = fh.read().splitlines()
    if not lines:
        continue
    if header is None:
        header = lines[0]
    rows += lines[1:]

# re-sort by druggability (col index 2) desc across all regions
def key(ln):
    try:
        return -float(ln.split("\t")[2])
    except Exception:
        return 0.0
rows.sort(key=key)

with open(out_path, "w") as out:
    out.write(header + "\n")
    out.write("\n".join(rows) + "\n")
print(f"[combine] wrote {out_path} ({len(rows)} pocket rows across {len(inputs)} regions)")
