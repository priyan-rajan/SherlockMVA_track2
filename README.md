# BUB1B N1002K — unified Nextflow pipeline (structure/druggability + MD)

One reproducible Nextflow workflow that runs **both** analyses built earlier for the
BUB1B `p.Asn1002Lys` variant:

- **Stage 1 — structural substrate & druggability** (env `bubdock`): fetch AlphaFold
  model of BUB1B (O60566), extract the pseudokinase region (750–1050), mutate Asn→Lys at
  1002, run fpocket druggability on full / WT-kinase / mutant-kinase, and a 97-pocket
  full-length survey.
- **Stage 2 — WT vs N1002K molecular dynamics** (env `bubmd`): prepare both systems,
  run 10 ns each (OpenMM, OBC2 implicit solvent, HMR/4 fs), then produce the
  pseudokinase packing/exposure report and the stability verdict.

The variant is `NM_001211.6:c.3006T>G p.Asn1002Lys` (15:40220612 T>G, GRCh38), in trans
with a truncating `p.Leu737Ter` allele. See `../../2026-10-23_bub1b-virtual-docking/` and
`../../2026-10-23_bub1b-md-stability/` for the motivating analyses.

## Quick start

```bash
cd ~/data/analysis/2026-10-24_bub1b-pipeline-nextflow
bash run_pipeline.sh -profile test      # fast smoke test (~minutes, tiny MD)
bash run_pipeline.sh                    # full run: 10 ns MD x2 (~12 h, 32 threads each)
bash run_pipeline.sh -profile structure_only   # Stage 1 only
bash run_pipeline.sh -resume            # resume after interruption
```

`run_pipeline.sh` auto-creates the two micromamba envs if missing, then runs
`nextflow run main.nf`.

## Outputs (`results/`)

| dir | contents |
|---|---|
| `01_structure/` | AlphaFold model + pLDDT |
| `02_regions/` | extracted WT kinase + N1002K mutant PDB |
| `03_druggability/` | `pocket_summary_*.tsv`, `fpocket_all_regions.tsv` |
| `04_pocket_survey/` | `pockets_annotated.tsv` (97-pocket full-length survey) |
| `05_md_prep/` | prepared solvated/parameterised systems |
| `06_md/` | `kinase_{WT,N1002K}.{dcd,final.pdb,log.csv,chk}` |
| `07_md_analysis/` | `pseudokinase_*`, `summary_stability.tsv`, `verdict.txt`, plots |
| `pipeline_trace.tsv`, `pipeline_report.html`, `run.log` | run metadata |

## Key parameters (`nextflow.config`)

| param | default | meaning |
|---|---|---|
| `acc` | `O60566` | UniProt accession |
| `kstart`/`kend` | 750 / 1050 | pseudokinase construct bounds |
| `mut_resid`/`mut_aa` | 1002 / LYS | variant |
| `env_struct`/`env_md` | `bubdock`/`bubmd` | micromamba envs |
| `run_md` | true | run Stage 2 |
| `md_ns` | 10 | production ns per state |
| `md_threads` | 32 | threads per MD run (2 run concurrently) |
| `save_ps` | 100 | trajectory frame spacing |

Override on the CLI, e.g. `bash run_pipeline.sh --md_ns 5 --md_threads 16`.

## Execution model

- Executor: `local` (no Docker/Singularity on this host).
- Tools invoked via `micromamba run -n <env> ...`.
- Stage 2 is gated: `ANALYZE` waits for all `MD` trajectories; the two systems run
  concurrently at `md_threads` each (needs ~2× threads ≤ node cores).
- Everything is checkpointed (`.chk`); use `-resume` to continue.

## Notes / caveats

- Only residues 750–1050 are modelled; the KARD domain (<766) is **not** included, so
  KARD/PP2A-B56 interactions cannot be assessed directly.
- MD is a single trajectory per state; treat the verdict as a structural hypothesis, not
  a ΔΔG. The decisive test remains wet-lab (BubR1 abundance/stability in patient cells).
- The mutant has a truncating co-allele; a missense-stabiliser cannot compensate for the
  null allele.

## Provenance
- Created: 2026-10-24
- `PI_SESSION_ID`: `01a16870-16e2-7157-82f7-9e3385823bb4`
- Reuses validated scripts from `../2026-10-23_bub1b-virtual-docking/bin/` and
  `../2026-10-23_bub1b-md-stability/bin/` (copied into `bin/`).
- Inputs: AlphaFold DB **AF-O60566-F1 v6** downloaded at runtime by `fetch_af.py`
  (UniProt O60566, BUB1B_HUMAN); variant `NM_001211.6:c.3006T>G p.Asn1002Lys`.
- Validation: Stage 1 executed for real (see `results/01..04`); full DAG incl. MD/Analyze
  validated with `nextflow run main.nf -stub-run` (all 13 processes COMPLETED).
  NOTE: `results/06_md` and `results/07_md_analysis` currently contain *empty placeholder*
  files from that stub run (see `STUB_PLACEHOLDERS.txt`); run the pipeline for real to fill them.
