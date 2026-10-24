#!/usr/bin/env python3
"""Extract a residue range [start,end] from a PDB file (all atoms, single chain)."""
import sys
from Bio.PDB import PDBParser, PDBIO, Select

pdb_in, pdb_out, start, end = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

class Region(Select):
    def accept_residue(self, res):
        _, num, _ = res.get_id()
        return start <= num <= end

p = PDBParser(QUIET=True)
s = p.get_structure("x", pdb_in)
io = PDBIO()
io.set_structure(s)
io.save(pdb_out, Region())
print(f"[extract_region] {pdb_in} -> {pdb_out} residues {start}-{end}")
