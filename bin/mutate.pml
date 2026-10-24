# PyMOL headless mutagenesis: WT -> mutant at one residue, with local repack.
# Invoked as: pymol -cq mutate.pml -- <in.pdb> <out.pdb> <resnum> <target_3letter>
import sys
from pymol import cmd

argv = sys.argv[1:]
if "--" in argv:
    argv = argv[argv.index("--") + 1:]
inp, outp, resnum, target = argv[0], argv[1], argv[2], argv[3]

cmd.load(inp, "prot")
cmd.wizard("mutagenesis")
cmd.get_wizard().do_select(f"{resnum}/")
cmd.get_wizard().set_mode(target)
cmd.get_wizard().apply()
cmd.set_wizard()
cmd.save(outp, "prot")
print(f"[mutate] {inp} residue {resnum} -> {target} written to {outp}")
