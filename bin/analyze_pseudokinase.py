#!/usr/bin/env python3
"""Variant-focused MD analysis for BUB1B N1002K: is the *pseudokinase domain*, and
specifically the 990-1015 surface that carries the variant and the PLK1 site T1008,
perturbed relative to WT?

Why this and not just global RMSD: the mechanistic claim we need to support is
"the N1002K substitution perturbs the C-terminal pseudokinase domain (UniProt
766-1050), the domain required for KARD phosphorylation -> PP2A-B56 recruitment ->
checkpoint silencing" (Cell Rep 2020, PMID 33207204). Global backbone RMSD cannot
speak to that. This script measures, WT vs N1002K, on the SAME trajectory set:

  1. domain-superposed global RMSD (766-1050) and local RMSD of segment 990-1015,
  2. per-residue RMSF grouped into functional sub-regions
     (750-765 juxtamembrane | 766-880 N-lobe | 881-989 | 990-1015 variant+PLK1 | 1016-1050),
  3. per-residue sidechain packing: number of heavy atoms within 4.5 A of each
     residue, and the WT->N1002K delta, i.e. a "packing-loss fingerprint",
  4. sidechain SASA of residue 1002 over time (burial vs exposure),
  5. per-residue helix occupancy (DSSP) around 1002.

Model scope: the construct is residues 750-1050, i.e. the pseudokinase domain plus
16 juxtamembrane residues. The KARD itself (<766) is NOT modelled -- state this as a
limitation; this script can only report on the domain that carries the variant.

Usage:
  python analyze_pseudokinase.py --dir results/02_md --out results/03_analysis [--ns-max 6]
"""
import argparse
import glob
import os

import matplotlib
import mdtraj as md
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True, help="dir with <label>.dcd and <label>.final.pdb")
ap.add_argument("--out", required=True)
ap.add_argument("--target", type=int, default=1002)
ap.add_argument("--domain", default="766-1050", help="superposition domain")
ap.add_argument("--window", default="990-1015", help="variant neighbourhood")
ap.add_argument("--contact-cut", type=float, default=0.45, help="nm, heavy-atom contact")
ap.add_argument("--save-ps", type=float, default=100.0, help="ps between saved frames (md_run.py --save-ps)")
ap.add_argument("--topdir", default=None,
               help="dir with <label>.pdb topologies (default: use --dir, then fall back to ../01_prep)")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

T = a.target
DLO, DHI = (int(x) for x in a.domain.split("-"))
WLO, WHI = (int(x) for x in a.window.split("-"))
REGIONS = [("juxta_750_765", 750, 765), ("Nlobe_766_880", 766, 880),
           ("core_881_989", 881, 989), ("variant_990_1015", 990, 1015),
           ("Cterm_1016_1050", 1016, 1050)]


def find_top(label, dirpath, topdir):
    cands = []
    if topdir:
        cands.append(os.path.join(topdir, label + ".pdb"))
    cands += [os.path.join(dirpath, label + ".final.pdb"),
              os.path.join(dirpath, label + ".pdb"),
              os.path.join(dirpath, "..", "01_prep", label + ".pdb")]
    for c in cands:
        if os.path.exists(c):
            return c
    return None


def load(dirpath, topdir):
    out = {}
    for dcd in sorted(glob.glob(os.path.join(dirpath, "*.dcd"))):
        label = os.path.basename(dcd)[:-4]
        top = find_top(label, dirpath, topdir)
        if top is None:
            print(f"[pk] skip {label}: no topology found")
            return out, label
        traj = md.load(dcd, top=top)
        if traj.n_frames < 3:
            print(f"[pk] skip {label}: {traj.n_frames} frames")
            continue
        out[label] = traj
    return out, None


trajs, missing = load(a.dir, a.topdir)
if not trajs:
    raise SystemExit(f"[pk] nothing to analyse in {a.dir} (missing topology for {missing})")

# common length so WT/mutant are compared over identical simulation time
n = min(t.n_frames for t in trajs.values())
for k in trajs:
    trajs[k] = trajs[k][:n]

stats, per_res, rmsf_tab, sasa_tab, ss_tab = {}, {}, [], [], []

for label, traj in trajs.items():
    ca_all = traj.topology.select("name CA")
    dom = np.array([i for i, at in enumerate(traj.topology.atoms)
                    if at.name == "CA" and DLO <= at.residue.resSeq <= DHI])
    traj.superpose(traj, frame=0, atom_indices=dom)

    resseq = np.array([traj.topology.atom(i).residue.resSeq for i in dom])
    rmsd = md.rmsd(traj, traj, 0, atom_indices=dom) * 10.0          # A

    ca = traj.atom_slice(ca_all)
    rmsf = md.rmsf(ca, ca, 0) * 10.0
    rseq_ca = np.array([r.resSeq for r in ca.topology.residues])

    # local RMSD of the variant window after domain superposition
    wsel = np.array([i for i, at in enumerate(traj.topology.atoms)
                     if at.name == "CA" and WLO <= at.residue.resSeq <= WHI])
    local_rmsd = md.rmsd(traj, traj, 0, atom_indices=wsel) * 10.0

    # heavy-atom packing: mean number of atoms within cut of each residue
    heavy = traj.topology.select("not element H")
    pairs = traj.topology.select_pairs(heavy, heavy)
    # keep only inter-residue pairs
    ri = np.array([traj.topology.atom(p[0]).residue.index for p in pairs])
    rj = np.array([traj.topology.atom(p[1]).residue.index for p in pairs])
    keep = ri != rj
    pairs, ri, rj = pairs[keep], ri[keep], rj[keep]
    res_of_index = [r.resSeq for r in traj.topology.residues]
    deg = np.zeros(traj.n_frames)
    # vectorised contact degree per residue, per frame
    deg_res = np.zeros((traj.n_frames, traj.topology.n_residues))
    for f in range(traj.n_frames):
        d = np.linalg.norm(traj.xyz[f][pairs[:, 0]] - traj.xyz[f][pairs[:, 1]], axis=1)
        near = d < a.contact_cut
        np.add.at(deg_res[f], ri[near], 1)
        np.add.at(deg_res[f], rj[near], 1)
    mean_deg = deg_res.mean(axis=0)
    sd_deg = deg_res.std(axis=0)

    # sidechain SASA of the variant residue
    # SASA in the context of the folded protein (whole-trajectory, per-residue)
    sasa_all = md.shrake_rupley(traj, mode="residue") * 100  # A^2, (frames, n_res)
    sasa = sasa_all[:, [r.index for r in traj.topology.residues if r.resSeq == T][0]]

    # helix occupancy per residue
    dssp = md.compute_dssp(traj, simplified=True)
    helix = (dssp == "H").mean(axis=0)

    ridx = {v: i for i, v in enumerate(res_of_index)}
    stats[label] = dict(nframes=n, ns=n * a.save_ps / 1000.0,
                        rmsd_dom=float(rmsd.mean()), rmsd_dom_sd=float(rmsd.std()),
                        rmsd_local=float(local_rmsd.mean()), rmsd_local_sd=float(local_rmsd.std()),
                        sasa_mean=float(sasa.mean()), sasa_last=float(sasa[-1]),
                        deg_target=float(mean_deg[ridx[T]]), deg_sd=float(sd_deg[ridx[T]]),
                        helix_1002=float(helix[ridx[T]]),
                        helix_window=float(helix[ridx[WLO]:ridx[WHI] + 1].mean()),
                        rmsf_dom=float(rmsf[(rseq_ca >= DLO) & (rseq_ca <= DHI)].mean()),
                        rmsf_target=float(rmsf[np.isclose(rseq_ca, T)][0]),
                        rmsf_window=float(rmsf[(rseq_ca >= WLO) & (rseq_ca <= WHI)].mean()),
                        regions={nm: float(rmsf[(rseq_ca >= lo) & (rseq_ca <= hi)].mean())
                                 for nm, lo, hi in REGIONS})
    per_res[label] = dict(deg=mean_deg, deg_sd=sd_deg, resseq=np.array(res_of_index), ss=helix)
    for r, v in zip(rseq_ca, rmsf):
        rmsf_tab.append((label, r, v))
    for f, v in enumerate(sasa):
        sasa_tab.append((label, f, v))

wt = next((k for k in trajs if "WT" in k and "N1002" not in k), None)
mu = next((k for k in trajs if "N1002" in k), None)

# ---------- outputs ----------
with open(os.path.join(a.out, "pseudokinase_summary.tsv"), "w") as fh:
    cols = ["system", "ns_analysed", "nframes", "rmsd_domain_A", "rmsd_domain_sd_A",
            "rmsd_variant_window_A", "rmsf_domain_A", "rmsf_1002_A", "rmsf_window_A",
            "packing_deg_1002", "packing_deg_sd", "sasa_1002_A2", "helix_occ_1002",
            "helix_occ_window"] + [r[0] for r in REGIONS]
    fh.write("\t".join(cols) + "\n")
    for k, s in stats.items():
        fh.write("\t".join([k, f"{s['ns']:.2f}", str(s["nframes"]),
                            f"{s['rmsd_dom']:.2f}", f"{s['rmsd_dom_sd']:.2f}",
                            f"{s['rmsd_local']:.2f}", f"{s['rmsf_dom']:.2f}",
                            f"{s['rmsf_target']:.2f}", f"{s['rmsf_window']:.2f}",
                            f"{s['deg_target']:.1f}", f"{s['deg_sd']:.1f}",
                            f"{s['sasa_mean']:.1f}", f"{s['helix_1002']:.2f}",
                            f"{s['helix_window']:.2f}"]
                           + [f"{s['regions'][r[0]]:.2f}" for r in REGIONS]) + "\n")

with open(os.path.join(a.out, "rmsf_by_region.tsv"), "w") as fh:
    fh.write("system\tregion\tresidues\trmsf_A\n")
    for k, s in stats.items():
        for nm, lo, hi in REGIONS:
            fh.write(f"{k}\t{nm}\t{lo}-{hi}\t{s['regions'][nm]:.2f}\n")

with open(os.path.join(a.out, "per_residue_rmsf.tsv"), "w") as fh:
    fh.write("system\tresSeq\trmsf_A\n")
    for row in rmsf_tab:
        fh.write(f"{row[0]}\t{row[1]}\t{row[2]:.3f}\n")

# packing-loss fingerprint
delta_rows = []
if wt and mu:
    dw, dm = per_res[wt], per_res[mu]
    delta = dm["deg"] - dw["deg"]
    for r, d in zip(dw["resseq"], delta):
        delta_rows.append((r, d))
    with open(os.path.join(a.out, "packing_delta_per_residue.tsv"), "w") as fh:
        fh.write("resSeq\tpacking_deg_WT\tpacking_deg_N1002K\tdelta\n")
        for r, d in delta_rows:
            i = int(np.where(dw["resseq"] == r)[0][0])
            fh.write(f"{r}\t{dw['deg'][i]:.1f}\t{dm['deg'][i]:.1f}\t{d:+.1f}\n")

# ---------- plots ----------
fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=False)
for k, s in per_res.items():
    ax[0].plot(s["resseq"], s["deg"], lw=1.0, label=k)
ax[0].axvline(T, color="k", ls="--", lw=0.8)
ax[0].set_ylabel("heavy-atom contacts / residue")
ax[0].set_title("Sidechain packing (4.5 Å), domain-superposed")
ax[0].legend(fontsize=8)

for k, s in per_res.items():
    ax[1].plot(s["resseq"], s["ss"], lw=1.0, label=k)
ax[1].axvline(T, color="k", ls="--", lw=0.8)
ax[1].set_ylabel("helix occupancy"); ax[1].set_xlabel("residue")
ax[1].set_title("DSSP helix occupancy")

if wt and mu:
    ax[2].bar([r for r, _ in delta_rows if WLO - 10 <= r <= WHI + 10],
              [d for r, d in delta_rows if WLO - 10 <= r <= WHI + 10])
    ax[2].axvline(T, color="k", ls="--", lw=0.8)
    ax[2].set_ylabel("Δcontacts (N1002K − WT)")
    ax[2].set_xlabel("residue"); ax[2].set_title("Packing-loss fingerprint around variant")
fig.tight_layout(); fig.savefig(os.path.join(a.out, "pseudokinase_packing.png"), dpi=140)

fig, ax = plt.subplots(figsize=(5, 4))
for k, s in stats.items():
    ax.bar(k, s["sasa_mean"])
ax.set_ylabel("residue 1002 sidechain SASA (Å²)")
ax.set_title("Variant sidechain burial (lower = more buried)")
fig.tight_layout(); fig.savefig(os.path.join(a.out, "sasa_1002.png"), dpi=140)

# ---------- verdict ----------
lines = []
if wt and mu:
    s_w, s_m = stats[wt], stats[mu]
    d_deg = s_m["deg_target"] - s_w["deg_target"]
    d_sasa = s_m["sasa_mean"] - s_w["sasa_mean"]
    d_rmsf = s_m["rmsf_window"] - s_w["rmsf_window"]
    d_loc = s_m["rmsd_local"] - s_w["rmsd_local"]
    d_helix = s_m["helix_window"] - s_w["helix_window"]
    flags = []
    if d_deg <= -1.0:
        flags.append(f"residue-1002 packing LOSS ({d_deg:+.1f} contacts)")
    elif d_deg >= 1.0:
        flags.append(f"residue-1002 packing GAIN ({d_deg:+.1f} contacts)")
    if d_sasa >= 15:
        flags.append(f"1002 sidechain becomes EXPOSED (ΔSASA {d_sasa:+.0f} Å²)")
    if d_rmsf >= 0.5:
        flags.append(f"variant window more flexible (ΔRMSF {d_rmsf:+.2f} Å)")
    if d_loc >= 0.5:
        flags.append(f"variant window moves relative to domain (ΔRMSD {d_loc:+.2f} Å)")
    if d_helix <= -0.15:
        flags.append(f"helix occupancy drops around 1002 (Δ{d_helix:+.2f})")
    verdict = ("; ".join(flags) if flags else
               "No local perturbation detected within this sampling — report as NEGATIVE "
               "for the pseudokinase-destabilisation hypothesis (do not over-claim)")
    lines = [f"WT vs N1002K, {s_w['ns']:.2f} ns analysed per system ({s_w['nframes']} frames)", "",
             f"domain backbone RMSD      WT {s_w['rmsd_dom']:.2f} Å | N1002K {s_m['rmsd_dom']:.2f} Å",
             f"variant-window RMSD       WT {s_w['rmsd_local']:.2f} Å | N1002K {s_m['rmsd_local']:.2f} Å",
             f"RMSF 990-1015             WT {s_w['rmsf_window']:.2f} Å | N1002K {s_m['rmsf_window']:.2f} Å",
             f"packing contacts @1002    WT {s_w['deg_target']:.1f} | N1002K {s_m['deg_target']:.1f}",
             f"SASA @1002 (sidechain)    WT {s_w['sasa_mean']:.0f} Å² | N1002K {s_m['sasa_mean']:.0f} Å²",
             f"helix occupancy 990-1015  WT {s_w["helix_window"]:.2f} | N1002K {s_m["helix_window"]:.2f}",
             "", f"VERDICT: {verdict}"]
else:
    verdict = "single-system run — no comparison possible"
    lines = [verdict]
with open(os.path.join(a.out, "pseudokinase_verdict.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"[pk] outputs in {a.out}")
