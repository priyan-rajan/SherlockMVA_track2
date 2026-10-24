#!/usr/bin/env python3
"""Implicit-solvent MD (Amber14 + OBC2) of a prepared protein with OpenMM (CPU).
Minimise -> Langevin equilibration -> production, writing DCD + CSV + checkpoint
so runs are restartable. Designed for CPU nodes (Threads platform property)."""
import argparse, os, sys, time
import openmm as mm
import openmm.app as app
import openmm.unit as u

ap = argparse.ArgumentParser()
ap.add_argument("--pdb", required=True)
ap.add_argument("--outprefix", required=True)
ap.add_argument("--ns", type=float, default=10.0)
ap.add_argument("--temp", type=float, default=300.0)
ap.add_argument("--threads", type=int, default=30)
ap.add_argument("--equil-ns", type=float, default=0.5)
ap.add_argument("--save-ps", type=float, default=100.0)
ap.add_argument("--resume", action="store_true")
ap.add_argument("--hmr", action="store_true",
                help="hydrogen mass repartitioning (enables 4 fs timestep)")
a = ap.parse_args()

def apply_hmr(system, topology, hmass=3.024):
    """Redistribute mass from heavy atoms onto their bonded H so a 4 fs timestep is
    stable with HBonds constraints (standard HMR). Conserves total mass."""
    n = 0
    for bond in topology.bonds():
        a1, a2 = bond
        e1 = a1.element.symbol if a1.element is not None else None
        e2 = a2.element.symbol if a2.element is not None else None
        if e1 == "H":
            h, heavy = a1.index, a2.index
        elif e2 == "H":
            h, heavy = a2.index, a1.index
        else:
            continue
        old = system.getParticleMass(h).value_in_unit(u.dalton)
        d = hmass - old
        system.setParticleMass(h, hmass * u.dalton)
        system.setParticleMass(
            heavy, (system.getParticleMass(heavy).value_in_unit(u.dalton) - d) * u.dalton)
        n += 1
    print(f"[md] HMR applied to {n} H atoms")

dt = 0.004 if a.hmr else 0.002  # ps; 4 fs with HMR + HBonds constraints
pdb = app.PDBFile(a.pdb)
ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")   # gbn2 = OBC2
system = ff.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic,
                         nonbondedCutoff=1.0*u.nanometer, constraints=app.HBonds,
                         soluteDielectric=1.0, solventDielectric=78.5)
if a.hmr:
    apply_hmr(system, pdb.topology)
integrator = mm.LangevinIntegrator(a.temp*u.kelvin, 1.0/u.picosecond, dt*u.picoseconds)
platform = mm.Platform.getPlatformByName("CPU")
props = {"Threads": str(a.threads)}
sim = app.Simulation(pdb.topology, system, integrator, platform, props)
sim.context.setPositions(pdb.positions)

os.makedirs(os.path.dirname(a.outprefix) or ".", exist_ok=True)
dcd = a.outprefix + ".dcd"
csv = a.outprefix + ".log.csv"
chk = a.outprefix + ".chk"
final = a.outprefix + ".final.pdb"

save_interval = max(1, int(round(a.save_ps / dt)))     # steps between saved frames
report_interval = 1000

if a.resume and os.path.exists(chk):
    sim.loadCheckpoint(chk); mode = "a"; print("[md] resumed from checkpoint")
else:
    print("[md] minimising..."); t0 = time.time()
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(a.temp*u.kelvin)
    mode = "w"
    print(f"[md] minimised in {time.time()-t0:.1f}s")

sim.reporters.append(app.DCDReporter(dcd, save_interval, append=(mode == "a")))
sim.reporters.append(app.StateDataReporter(csv, report_interval, append=(mode == "a"),
    step=True, time=True, potentialEnergy=True, temperature=True, speed=True))
sim.reporters.append(app.CheckpointReporter(chk, save_interval))

n_equil = int(round(a.equil_ns * 1000.0 / dt))   # ns->ps / (ps/step)
n_prod = int(round(a.ns * 1000.0 / dt))
print(f"[md] equil {a.equil_ns} ns ({n_equil} steps) + production {a.ns} ns ({n_prod} steps), dt={dt} ps, threads={a.threads}")
t0 = time.time()
total = n_equil + n_prod
done = 0
chunk = 25000
while done < total:
    n = min(chunk, total - done)
    sim.step(n)
    done += n
    el = time.time() - t0
    print(f"[md] {done}/{total} steps  {done/total*100:4.1f}%  elapsed {el/60:6.1f} min", flush=True)

with open(final, "w") as fh:
    app.PDBFile.writeFile(sim.topology, sim.context.getState(getPositions=True).getPositions(),
                          fh, keepIds=True)
print(f"[md] done {a.outprefix} in {(time.time()-t0)/60:.1f} min")
