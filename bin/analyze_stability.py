#!/usr/bin/env python3
"""Compare WT vs N1002K kinase-domain MD. Quantifies global and *local* stability
around residue 1002: backbone RMSD, per-residue RMSF, Cα contact count and DSSP
secondary-structure content. Writes a summary TSV, per-residue TSVs and plots, and
prints a plain-language verdict on whether N1002K locally destabilises the domain."""
import argparse, glob, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True, help="dir with <label>.dcd + <label>.final.pdb")
ap.add_argument("--out", required=True)
ap.add_argument("--target", type=int, default=1002)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
T = a.target

systems = {}
for dcd in sorted(glob.glob(os.path.join(a.dir, "*.dcd"))):
    label = os.path.basename(dcd)[:-4]
    top = os.path.join(a.dir, label + ".final.pdb")
    if not os.path.exists(top):
        print(f"[analyze] skip {label}: no topology"); continue
    traj = md.load(dcd, top=top)
    if traj.n_frames < 2:
        print(f"[analyze] skip {label}: only {traj.n_frames} frame(s)"); continue
    ca = traj.atom_slice(traj.topology.select("name CA"))
    rmsd = md.rmsd(ca, ca, 0) * 10.0          # nm -> Å
    rmsf = md.rmsf(ca, ca, 0) * 10.0          # nm -> Å
    resseq = np.array([r.resSeq for r in ca.topology.residues])

    def rmsf_rng(lo, hi):
        m = (resseq >= lo) & (resseq <= hi)
        return float(np.mean(rmsf[m])) if m.any() else float("nan")

    ti = np.where(resseq == T)[0]
    if len(ti):
        ti = ti[0]
        xyz = ca.xyz                              # nm, (frames, nca, 3)
        d = np.linalg.norm(xyz - xyz[:, ti:ti + 1, :], axis=2)
        contacts = (d < 0.8).sum(axis=1) - 1      # Cα neighbours within 8 A
        mean_contacts = float(np.mean(contacts))
    else:
        contacts = None; mean_contacts = float("nan")

    dssp = md.compute_dssp(traj, simplified=True)  # H/E/C
    systems[label] = dict(
        label=label, time=traj.time / 1000.0, rmsd=rmsd, rmsf=rmsf, resseq=resseq,
        mean_rmsd=float(np.mean(rmsd)), final_rmsd=float(rmsd[-1]), max_rmsd=float(np.max(rmsd)),
        mean_rmsf=float(np.mean(rmsf)),
        rmsf_target=rmsf_rng(T - 3, T + 3), rmsf_far=rmsf_rng(T + 40, T + 60),
        mean_contacts=mean_contacts, contacts=contacts,
        frac_helix=float(np.mean(dssp == "H")), frac_sheet=float(np.mean(dssp == "E")),
        nframes=traj.n_frames)

# --- write outputs ---
with open(os.path.join(a.out, "summary_stability.tsv"), "w") as fh:
    fh.write("system\tnframes\tmean_rmsd_A\tfinal_rmsd_A\tmax_rmsd_A\tmean_rmsf_A\t"
             "rmsf_1002_pm3_A\trmsf_distal_A\tmean_CA_contacts_1002\tfrac_helix\tfrac_sheet\n")
    for s in systems.values():
        fh.write(f"{s['label']}\t{s['nframes']}\t{s['mean_rmsd']:.2f}\t{s['final_rmsd']:.2f}\t"
                 f"{s['max_rmsd']:.2f}\t{s['mean_rmsf']:.2f}\t{s['rmsf_target']:.2f}\t{s['rmsf_far']:.2f}\t"
                 f"{s['mean_contacts']:.1f}\t{s['frac_helix']:.3f}\t{s['frac_sheet']:.3f}\n")
for s in systems.values():
    with open(os.path.join(a.out, f"rmsf_{s['label']}.tsv"), "w") as fh:
        fh.write("resSeq\trmsf_A\n")
        for r, v in zip(s["resseq"], s["rmsf"]):
            fh.write(f"{r}\t{v:.3f}\n")

# --- plots ---
fig, ax = plt.subplots(figsize=(7, 4))
for s in systems.values():
    ax.plot(s["time"], s["rmsd"], label=s["label"], lw=1.2)
ax.set_xlabel("time (ns)"); ax.set_ylabel("backbone RMSD (Å)"); ax.set_title("Kinase domain RMSD")
ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(a.out, "rmsd_vs_time.png"), dpi=140)

fig, ax = plt.subplots(figsize=(9, 4))
for s in systems.values():
    ax.plot(s["resseq"], s["rmsf"], label=s["label"], lw=1.0)
ax.axvline(T, color="k", ls="--", lw=0.8); ax.text(T, ax.get_ylim()[1]*0.95, " N1002K", fontsize=8)
ax.set_xlabel("residue"); ax.set_ylabel("per-residue RMSF (Å)"); ax.set_title("Cα RMSF (WT vs N1002K)")
ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(a.out, "rmsf_vs_residue.png"), dpi=140)

# --- verdict ---
wt = next((s for k, s in systems.items() if "WT" in k and "N1002" not in k), None)
mu = next((s for k, s in systems.items() if "N1002" in k), None)
verdict = "insufficient data for WT/mutant comparison"
if wt and mu:
    rf_ratio = mu["rmsf_target"] / wt["rmsf_target"] if wt["rmsf_target"] else float("nan")
    dc = mu["mean_contacts"] - wt["mean_contacts"]
    drmsd = mu["final_rmsd"] - wt["final_rmsd"]
    flags = []
    if rf_ratio == rf_ratio and rf_ratio > 1.5: flags.append("local flexibility ↑>1.5× at 1002")
    if dc < -1.0: flags.append("loss of 1002 packing contacts")
    if drmsd > 2.0: flags.append("mutant backbone diverges (ΔRMSD>2 Å)")
    verdict = ("N1002K shows LOCAL DESTABILISATION: " + "; ".join(flags)) if flags else \
              (f"No clear local destabilisation (rmsf ratio {rf_ratio:.2f}, Δcontacts {dc:+.1f}, "
               f"ΔRMSD {drmsd:+.2f} Å)")
with open(os.path.join(a.out, "verdict.txt"), "w") as fh:
    fh.write(verdict + "\n")
print("[analyze] " + verdict)
print(f"[analyze] outputs in {a.out}")
