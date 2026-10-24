#!/usr/bin/env python3
"""Prepare a raw (AlphaFold / in-silico-mutated) PDB fragment for MD.
Fix missing atoms/residues, strip heterogens, protonate at pH 7.0. All-atom, no solvent
(implicit-solvent MD)."""
import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile

ap = argparse.ArgumentParser()
ap.add_argument("--inp", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--ph", type=float, default=7.0)
a = ap.parse_args()

fixer = PDBFixer(filename=a.inp)
fixer.removeHeterogens(keepWater=False)
fixer.findMissingResidues()
fixer.missingResidues = {}          # fragment: do NOT invent disordered loops
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(a.ph)

with open(a.out, "w") as fh:
    PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

n_atoms = sum(1 for _ in open(a.out) if _.startswith(("ATOM", "HETATM")))
resids = sorted({int(l[22:26]) for l in open(a.out) if l.startswith("ATOM")})
print(f"[prep] {a.out}: {n_atoms} atoms, residues {resids[0]}-{resids[-1]} ({len(resids)} res)")
