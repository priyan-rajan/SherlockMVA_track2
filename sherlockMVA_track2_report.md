# Track 2 — Drug Repurposing Report

## Rescuing Residual BubR1 Function in Candidate Compound-Heterozygous Mosaic Variegated Aneuploidy

*A mechanism-driven repurposing strategy for a child carrying `c.2210T>G (p.Leu737Ter)` and
`c.3006T>G (p.Asn1002Lys)` in BUB1B*

> **Status:** hypothesis for further investigation — **not** evidence that any medicine works.
> Structure: the disease (§1) → the two variants (§2) → molecular dynamics (§3) → docking
> (§4) → proposed therapies (§5) → diagnostics (§6) → acknowledgements (§7) → references (§8).
> Every external claim is traceable to a PMID (see `evidence/` and
> `analysis/2026-10-24_track2-reconciliation/results/verification.txt`).
> Full citations are listed in **§8 References**; PMIDs are also given inline at each claim.

---

## 1. Mosaic Variegated Aneuploidy

Mosaic Variegated Aneuploidy (MVA1, OMIM 257300) is caused by failure of the **spindle
assembly checkpoint (SAC)** — the surveillance system that holds mitosis open until every
chromosome is correctly attached to the mitotic spindle. The SAC's central scaffold is
**BubR1**, encoded by *BUB1B* on chromosome 15q15. The protein's checkpoint function is
largely independent of its kinase activity; its C-terminal domain is a degenerate
**pseudokinase** that instead acts as a protein–protein-interaction module, relaying KARD
phosphorylation to PP2A-B56 for checkpoint silencing and chromosome alignment (PMID 33207204).

When BubR1 is limiting, cells slip through metaphase before attachments are corrected.
Daughter cells inherit the wrong chromosome number — a state that is *variegated* (different
chromosome errors in different cells) and *mosaic* (present in a fraction of cells). The
clinical consequences are microcephaly, developmental delay and growth failure, together with
a sharply elevated childhood cancer risk (Wilms tumour, rhabdomyosarcoma, leukaemia).

**The therapeutic logic follows directly from the biology.** Because the defect is a
*shortage* of BubR1, the useful direction is to **raise the amount of functional BubR1** —
not to inhibit anything. Figure 1 sets out the whole argument on one page.

![Figure 1 — One genotype, three hypotheses, one definitive fix. Two different *BUB1B* alleles
converge on a single quantitative defect (too little BubR1, so the checkpoint cannot hold),
which opens three testable hypotheses — A readthrough, B proteasome inhibition, C
acetylation-stability — and one definitive editing route. (`figures/fig1_the_story.png`)](figures/fig1_the_story.png)

**Figure 1.** One genotype, three hypotheses, one definitive fix. Only **A** addresses the defect
itself; **C** cannot create protein, only keep what is made (§5.3); editing is the definitive
future fix but is not a repurposing candidate.

---

## 2. The Two Variants

| Allele | GRCh38 | cDNA (MANE `NM_001211.6` / `ENST00000287598`) | Protein | Location | Molecular consequence |
|---|---|---|---|---|---|
| **A** | 15:40209701 T>G | c.2210T>G | p.Leu737Ter | Exon 17/23 (last 3′ junction at c.2957) | **UGA (opal) stop** → NMD-competent transcript → little or no full-length protein |
| **B** | 15:40220612 T>G | c.3006T>G | p.Asn1002Lys | Exon 23/23, pseudokinase domain (766–1050) | Missense → protein produced; function partially retained (VUS) |

Both variants are heterozygous (`0/1`), PASS, with read depths 21/25 and 15/13 (total DP 46 and
28). They come from the Sentieon/GATK small-variant VCF for the singleton proband `WGS_EX2312012`
(GRCh38).

**Allele A** is ClinVar **Pathogenic/Likely_pathogenic** for MVA1 (multiple submitters, no
conflicts) and has a gnomAD v4.1 frequency of 7.85×10⁻⁵ (AC 115 / 1,461,846). SpliceAI max DS is
0.03, so it is not splice-altering. It is a well-established loss-of-function MVA allele.

**Allele B** is not in ClinVar (a VUS). Its supporting in-silico evidence is AlphaMissense 0.923,
ESM-1v +0.465, REVEL 0.472 (just below the 0.5 support threshold), PhyloP 3.691 and gnomAD
frequency 5×10⁻⁷ (1 allele in 1,461,878).

**Concordance with the Track 1 analysis.** The same pair was nominated independently by the
causal-variant workflow. Exomiser ranked *BUB1B* first under the dominant model and second under
the recessive model. No other MVA-panel gene (*CEP57*, *TRIP13*, *MAD1L1*, *SMC5*) carried a
plausible variant, and a genome-wide SV/CNV scan (Manta, CNVkit, TIDDIT) found no credible MVA-gene
structural or copy-number event. This does not resolve phasing (see below).

> **Phasing is not established.** **These are *candidate* compound heterozygous variants, not a
> phased genotype.** The call set does not determine whether c.2210T>G and c.3006T>G lie on
> opposite (*trans*) or the same (*cis*) homologues, and no parental segregation or long-read
> phasing result was available for this report. *Trans* is the expected configuration for a
> recessive MVA1 phenotype and the one this report works from, but it is an inference from the
> phenotype, **not a measurement** — and it matters: if both variants sit in *cis*, the other
> allele would be wild-type, the genotype would not explain the disease, and this analysis would
> need re-interpretation. **Confirmatory phasing (parental genotypes, or long-read sequencing
> across the locus) is the first thing to close.** Everything downstream is conditional on it.

**Allele A is the null arm.** The **premature termination codon (PTC)** — the opal (UGA) stop
that `c.2210T>G` puts in place of leucine 737, cutting the protein off at 737 of its 1050 residues
instead of letting the ribosome read through to the normal stop — sits **747 nt upstream of the
last exon–exon junction**, well beyond the ~50–55 nt rule. A PTC in that position marks the
transcript for **nonsense-mediated decay (NMD)**, so the mRNA is degraded rather than translated:
the transcript is **NMD-competent** and contributes essentially no protein
(`evidence/ptc_context.md`, `evidence/ptc_context.tsv`).

**Allele B is the hypomorphic arm — and, under the candidate compound-heterozygous model, the
child's only source of BubR1.** The
missense sits inside the pseudokinase domain, six residues from the PLK1 site T1008. Known
MVA1 missense alleles cluster in this same domain (814, 844, 909, 921, 1012), all annotated
as compound heterozygous with a nonsense allele — matching the patient's *candidate*
architecture — which
makes N1002K a plausible hypomorph by analogy to five disease alleles. **The exact residual
amount is not yet measured**, and that single number converts this proposal into arithmetic
(§5.1, §5.6).

### 2.1 The design principle: do not inhibit BUB1B

The instinctive oncology move is to inhibit the protein, and there is a real rationale for it:
aneuploid cells already carry a broken checkpoint and an imbalanced proteome, so stripping what
checkpoint function they have left may kill them preferentially (§5.2 relies on the same idea).
The problem is that **BubR1 is not tumour-specific** — it is a housekeeping protein in every
dividing cell — so the same drug also lowers checkpoint fidelity in cells that are currently
normal, and each of their divisions becomes a chance to create a *new* aneuploid cell. The net
effect in a child is not "kill the aneuploid cells" but "kill some of them and make more", which
is why *BUB1B* inhibition is excluded as an anti-target rather than proposed as a therapy.

The same reasoning eliminates the entire anti-target set:

- **MPS1, PLK1 and Aurora B inhibitors** — weaken an already weak checkpoint.
- **Taxanes** — increase the very mis-segregation the checkpoint cannot catch.
- **PP2A activators** — accelerate checkpoint silencing, the process already running too fast.
- **CC-90009 and other cereblon modulators** — degrade the protein that is already limiting.
- **NMD inhibition alone (e.g. SMG1 inhibitor)** — releases the PTC transcript *without*
  readthrough, yielding only the 737-aa fragment (it retains the KARD and kinetochore-binding
  N-terminus), a potential **dominant-negative**; rejected on mechanism. NMD inhibition is
  acceptable only *alongside* readthrough, never instead of it.

---

## 3. Molecular Dynamics: Method and Analysis

The mutation's structural consequence has to be established before any strategy is designed around
it. MD was run to answer one question directly: **is N1002K actually destabilising the
pseudokinase domain, or is it a benign substitution?** Full detail:
`analysis/2026-10-23_bub1b-md-stability/`.

### 3.1 Method

| Item | Setting |
|---|---|
| Systems | pseudokinase region **750–1050** (301 residues, ~4.9k atoms), **WT vs N1002K** |
| Input | extracted from AlphaFold model **AF-O60566-F1 (v6)**; mutant built by PyMOL mutagenesis with local repack in the docking analysis upstream; protonated and fixed with **pdbfixer** |
| Engine | **OpenMM**; force field **amber14-all + OBC2 (GBSA) implicit solvent** |
| Parameters | 1.0 nm non-bonded cutoff, HBonds constraints, **4 fs timestep with HMR**, Langevin thermostat **300 K**, friction 1/ps |
| Protocol | energy minimisation (2000 iterations) → 0.5 ns equilibration → **10 ns production**; frames every 100 ps |
| Sampling achieved | **10.5 ns per system, 105 frames**; WT and mutant run concurrently at 32 threads each (measured throughput optimum) |
| Analysis | `mdtraj`: backbone RMSD, Cα RMSF, heavy-atom packing (other residues' non-H atoms within 4.5 Å, per residue), residue-1002 SASA, DSSP helix/sheet occupancy; all values are **means over the 105 frames of the 10.5 ns run** |

### 3.2 Analysis and result — local packing loss, not unfolding

| Quantity (WT → N1002K) | WT | N1002K | Δ |
|---|---|---|---|
| Domain backbone RMSD (766–1050) | 2.62 Å | 2.62 Å | 0.00 |
| Variant-window RMSD (990–1015) | 0.85 Å | 1.49 Å | **+0.64 Å** |
| Variant-window RMSF (990–1015) | 0.97 Å | 1.48 Å | **+0.52 Å** |
| Packing contacts at residue 1002 (mean, 10.5 ns) | 85.2 | 69.4 | **−15.8** |
| SASA at residue 1002 (mean, 10.5 ns) | 32 Å² | 112 Å² | **+80 Å²** |
| SASA at residue 1002 (frame range, 10.5 ns) | 11–62 Å² | 73–163 Å² | **no overlap** |
| Helix occupancy (990–1015) | 0.79 | 0.76 | −0.03 |

**Verdict: N1002K causes local structural destabilisation without global unfolding.** The
domain folds to the same global RMSD in both states and secondary structure is preserved, but
the variant window becomes measurably more mobile (+0.64 Å RMSD, +0.52 Å RMSF), the side chain
loses about sixteen packing contacts, and it goes from buried (32 Å²) to exposed (112 Å²) **on
average**.
The affected patch is contiguous and sits immediately N-terminal to the **PLK1 site T1008**.
The signal was already present at 2.7 ns and persisted at 10.5 ns.

> **These are trajectory means, not single numbers — and that matters for how the animation
> reads.** Residue-1002 SASA is measured every frame and moves a lot: WT runs **11–62 Å²**, N1002K
> runs **73–163 Å²** across the same 105 frames. The two ranges do not overlap, which is why the
> mean difference is trustworthy, but an individual frame can read 18 Å² or 78 Å². When the MD
> clip's live counter is quoted, it is quoting one frame of that range, not the 32 → 112 Å² mean.

![Figure 2 — Four-panel MD comparison, WT vs N1002K, all panels means over 10.5 ns: packing
contacts at residue 1002 (85 → 69), SASA at residue 1002 (32 → 112 Å²), RMSF of the variant
window 990–1015 (0.97 → 1.48 Å), and local RMSD of that window (0.85 → 1.49 Å).
(`figures/fig3_md_packing_4panel.png`)](figures/fig3_md_packing_4panel.png)

**Figure 2.** WT (blue) vs N1002K (red). The domain does not unfold (global backbone RMSD is
2.62 Å in both states, not shown); a local window around the variant loosens and becomes exposed.
Bars are **10.5 ns means**, which is why they have no error bars — the frame-to-frame SASA swing
is large (WT 11–62 Å², N1002K 73–163 Å², §3.2). Data:
`results/03_analysis/pseudokinase_summary.tsv`.

![Figure 3 — Per-residue packing change (WT minus N1002K) across the pseudokinase domain. The
loss is concentrated at the variant site rather than spread across the domain.
(`figures/fig3b_packing_per_residue.png`)](figures/fig3b_packing_per_residue.png)

**Figure 3.** Per-residue packing delta. Data: `packing_delta_per_residue.tsv`.

### 3.3 Interpretation and caveats

Two consequences follow. First, the mutant is *destabilised*. Second, the destabilisation is
*local and partial*, which is consistent with a **hypomorph** — exactly what allele B needs to be
for this genotype to be survivable at all.

**Caveats, stated plainly:** a single trajectory per state; implicit solvent; no KARD domain
in the model (so inter-domain effects are invisible); part of the ΔSASA is simply the larger
lysine side chain; and this is not a ΔΔG. Replicates and free-energy calculations, plus
wet-lab measurement of BubR1 abundance and stability, would be needed to make it decisive.

---

## 4. Docking: A Negative Result That Shapes the Strategy

MD answers "is the mutant destabilised?". The complementary question — and the one that
determines whether a small molecule could be designed at all — is **"is there anything to bind?"**
Pocket detection was run to find out.

### 4.1 Method

Full workflow and provenance: `analysis/2026-10-23_bub1b-virtual-docking/`.

- **Structure** — AlphaFold DB model **AF-O60566-F1 (v6)** for UniProt **O60566**
  (BUB1B_HUMAN, 1050 aa); per-residue pLDDT retrieved alongside.
- **Regions** — full-length protein, and the C-terminal region **750–1050** (pseudokinase
  domain plus juxtamembrane residues).
- **Mutant** — WT → **N1002K** by PyMOL mutagenesis with local repack
  (`pymol-open-source 3.1.0`); residue identity verified in the output.
- **Pocket detection** — **fpocket 4.2.3** on full-length, kinase WT and kinase N1002K;
  pockets parsed and ranked by druggability score, and flagged when they line residue 1002.
- **Tools** — openbabel 3.2.1, pdb-tools 2.7.0, biopython 1.88, Python 3.11.17; orchestrated
  by Nextflow 25.10.4 in the `bubdock` environment. 151 pocket rows combined across regions.

The question we asked was: **does the missense create a pocket that a small molecule could
occupy to stabilise the protein?** (A *ligand-docking* campaign against *BUB1B*
was not run, deliberately: docking a chaperone/proteasome/checkpoint library into BubR1 would
address the wrong target entirely — those compounds act on their own proteins.)

### 4.2 Result — no druggable pocket

- **The site itself is confident**, not a disordered artefact: Cα pLDDT at residue 1002 is
  **91.1** (residues 1000–1005: 93.8 / 91.6 / 91.1 / 88.0 / 80.4 / 82.2).
- **Nothing strongly druggable anywhere.** Best score across all 97 full-length pockets is
  **0.169** (pocket 9, 1514 Å³); the best kinase-region pocket is **0.132** — both low on the
  usual ≥0.5 druggability scale, consistent with a scaffold/PPI surface.
- **The mutation creates no new cavity.** WT and N1002K kinase regions share the *identical*
  top pocket (0.132). The pocket that lines residue 1002 scores **0.005** in both states and
  changes only marginally (WT pocket 5 vs mutant pocket 6, 374.3 vs 374.4 Å³).

![Figure 4 — fpocket druggability across BUB1B regions. The pocket lining the variant site
scores 0.005, far below the ~0.5 threshold usually taken as ligandable; WT and N1002K are
indistinguishable. (`figures/fig4_fpocket_no_pocket.png`)](figures/fig4_fpocket_no_pocket.png)

**Figure 4.** Druggability survey. Bars are the best-scoring pocket per region; the pocket that
lines residue 1002 is essentially undruggable and the mutation does not change it. Data:
`fpocket_all_regions.tsv` (151 pocket rows).

### 4.3 What this rules out — and what it redirects to

Taken with §3, the two structural results tell one consistent story: **N1002K loosens packing and
exposes a side chain, but leaves no cavity to fill.** The negative result therefore eliminates the
**"direct stabiliser" arm** — a chemical chaperone that binds N1002K and props it up is not a
realistic starting point — and it eliminates the framing in which docking is performed *against
BubR1* at all.

It redirects instead to the two levers that do not require binding the mutant: **how much
protein the cell makes** (§5.1) and **how the cell handles the aneuploid state** (§5.2).
A negative result that removes a whole strategy is worth reporting as one — it is the reason
this report does not propose a BubR1-binding drug.

---

## 5. Proposed Therapies

Three hypotheses — **A** readthrough (§5.1), **B** proteasome inhibition (§5.2), **C**
acetylation-stability (§5.3), ordered by **mechanism strength against risk in a child**. Two
further sections are deliberately *not* therapeutic proposals but investigational routes: a
delivery mode (§5.4) and gene editing (§5.5), which would be required to make any of this
permanent. §5.6 is the experiment that prices the whole programme. The diagnostic content for the
next affected family is in §6.
Figure 5 is the trade-off view.

![Figure 5 — Hypothesis trade-off: mechanism strength versus toxicity risk in a paediatric patient.
Readthrough sits in the high-benefit / manageable-risk corner; proteasome inhibition sits in
the high-risk corner and is therefore last resort. (`figures/fig5_tier_tradeoff.png`)](figures/fig5_tier_tradeoff.png)

**Figure 5.** Mechanism strength vs toxicity risk. The ordering follows this trade-off, not
enthusiasm for any single drug: **C** ranks below **A** because it cannot create protein, but it
is the least risky lever and multiplies A rather than competing with it.

The framework also generalises: the "one null arm plus one hypomorphic arm" architecture
(assumed here as a candidate, pending phasing — §2)
applies to other compound-heterozygous *BUB1B* patients, and the readthrough lever transfers to
nonsense variants in other MVA genes (*CEP57*, *TRIP13*, *MAD1L1*, *SMC5*), where the same
arithmetic — residual % versus the ~50% threshold — decides feasibility.

### 5.1 Hypothesis A — Codon readthrough: the mechanism fix

The `p.Leu737Ter` transcript carries an **opal (UGA)** codon, and readthrough there restores
full-length, NMD-resistant BubR1. This is the only lever upstream of *all three* defective arms
of the checkpoint (mitotic checkpoint complex, error correction, premature chromatid
separation).

**This is proposed as a drug *class*, not a single molecule.** Codon specificity, potency,
safety and development status differ sharply between agents, so the defensible move is to
**screen the class in patient-derived cells and let the data choose**. Ranking detail and
PMIDs: `docs/02_readthrough_modalities.md`.

| Candidate | Fit for **UGA** (our PTC) | Status | Key evidence |
|---|---|---|---|
| **ELX-02** | Strongest — PTC-selective, UGA-competent decoding-site agent | Clinical-stage (Phase 2 in nonsense CF/Alport); **not approved** | Readthrough *without* native-stop readthrough (PMID 32376628); modest, heterogeneous rescue across 206 patient organoids (PMID 42728176) |
| **Amlexanox** | Good — readthrough **plus** NMD inhibition (dual action) | **Approved** (topical, aphthous ulcers); oral form in metabolic trials | Repurposing candidate: approved, orally available, dual mechanism |
| **Ataluren** | **Weak — see caveat** | Conditionally approved (EU) for nonsense DMD; status under review | UAA-biased; in vivo rescued **only** the UAA mutant, with **no effect at UAG or UGA** (PMID 32343037) |
| **Gentamicin** | Reads through all three stops, incl. UGA | Approved antibiotic; decades of paediatric use | Paediatric human precedent: 5 children, IV, protein restored (PMID 35234826) |

**Ataluren caveat (important).** Readthrough is codon-dependent, and termination strength is
UAA > UAG > **UGA** (PMID 32343037). Our PTC is the hardest class, so ataluren is included as a
class member but **ranks below ELX-02 and amlexanox** for this genotype and should not be the
lead. We state this because a reader who checks the source will find it.

**Two non-drug modalities in the same class:**

- **Low-dose potentiator construct (the novel design).** eRF3/GSPT1 degraders do **not** read
  through alone — they *potentiate* aminoglycosides (PMID 42907672, 33764477), permitting lower
  and less toxic doses of each agent. A low-dose decoding-site agent combined with a low-dose
  release-factor potentiator is a design no current programme offers.
- **Engineered suppressor tRNA (scalability horizon).** Non-viral LNP delivery of chemically
  modified sup-tRNA now works in vivo (PMID 42658930) — a genotype-agnostic route if small
  molecules underdeliver.

**The arithmetic that makes this a proposal rather than a slogan.** In graded BubR1 knockdown
(PMID 18932004), premature chromatid separation appears at *every* level, but **aneuploidy
only below ~50% residual BubR1**. Therefore:

```
required readthrough  =  50%  −  (measured BubR1 in patient cells, % of control)
```

![Figure 6 — Why the target is arithmetic. Schematic aneuploidy probability as a function of
residual BubR1: the phenotype switches below ~50% (PMID 18932004), so the therapeutic target is
that threshold minus the patient's measured level. The patient's value is not yet measured.
(`figures/fig2_dosage_threshold.png`)](figures/fig2_dosage_threshold.png)

**Figure 6.** The dosage-threshold argument. The curve is schematic, drawn from the graded
knockdown result in PMID 18932004; the patient's residual BubR1 is the missing number, not a
measured value.

A few percent of full-length protein may be enough to cross the threshold. That is a
falsifiable target — unlike "restore BubR1".

**Safety framing.** The relevant question is the **selectivity index** (PTC readthrough ÷
native-stop readthrough), not toxicity in the abstract; reduced toxicity has been an explicit
design goal of the aminoglycoside series since 2009 (NB54 → NB124 → ELX-02). A
genotype-specific index — Δ(full-length BubR1) ÷ Δ(truncated 737-aa fragment) — should
accompany any readthrough dose, because NMD inhibition raises both products.

**Figure 7 is the codon-level view.** The same two codons that constrain readthrough also
define what a permanent genetic fix would have to do, so it is introduced here and re-used in
§5.5.

![Figure 7 — Two alleles, two different editing problems. Allele A: an adenine base editor
converts the opal stop TGA to TGG (Trp) with no double-strand break. Allele B: AAT→AAG is a
transversion requiring prime editing. (`figures/fig6_editing_codons.png`)](figures/fig6_editing_codons.png)

**Figure 7.** The two alleles at codon level: a *readthrough* problem here, an *editing* problem
in §5.5. The asymmetry — one base edit versus one transversion —
is what makes allele A the tractable target and allele B the harder prize.

### 5.2 Hypothesis B — Proteasome inhibitors: chemoprevention and last resort

**Least favoured; retained because the biology is real, not because the safety case is good.**
Aneuploid cells carry an imbalanced proteome and a chronic proteotoxic burden — aneuploid
embryos show large-scale transcriptome dysregulation (PMID 30297919), aneuploidy imposes a
protein-degradation and quality-control load (PMID 31068692), and specific degradation
dependencies (e.g. UBE2H) have been mapped (PMID 39247952, 42094535). The argument is
**synthetic lethality**: the aneuploid cell's coping mechanism becomes its vulnerability.

- **Bortezomib (proteasome inhibitor, FDA-approved).** Ultra-low, intermittent dosing would be
  needed to nudge turnover without systemic toxicity. **Hard caveats:** at chemotherapy doses
  it is neurotoxic and myelosuppressive; a proteasome inhibitor in a genome-instability child
  is a **cytotoxic** proposition; and the claim that proteasome inhibition *raises* mutant
  BubR1 is **a testable hypothesis, not established** — the MD data (§3) argue the mutant is
  destabilised and exposed, not stabilised. Defensible use: chemoprevention in
  aneuploid-burdened tissue, in patient-derived cells first.
- **Genotype caveat — the patient's LZTR1 p.Tyr748Ter may blunt proteasome-inhibitor efficacy.**
  The WGS also carries a **heterozygous LZTR1 truncating variant** (chr22:20996720 C>G,
  `c.2244C>G`, `p.Tyr748Ter`; LP, VAF ~48%; a secondary finding, Noonan-spectrum /
  LZTR1-related schwannomatosis). LZTR1 is the substrate **adaptor** of a CUL3-RING E3
  ubiquitin ligase that tags RAS-family GTPases (KRAS/NRAS/HRAS/MRAS) for **proteasomal
  degradation** — it supplies the *specificity*, upstream of the proteasome. A proteasome
  inhibitor blocks the **downstream executioner**, so it cannot restore LZTR1 function; on the
  RAS axis the two push the **same direction** — LZTR1 haploinsufficiency already impairs RAS
  turnover, and proteasome inhibition reduces turnover further, **raising RAS/MAPK output**.
  Because RAS/MAPK activation is a recognised, context-dependent mechanism of **resistance to
  proteasome inhibitors**, the LZTR1 variant is a plausible cause of **reduced bortezomib /
  carfilzomib efficacy** in this patient. This is a hypothesis, not established, and is a further
  reason to answer it empirically — RAS-GTP and pERK ± proteasome inhibitor in patient-derived
  cells — before any *in vivo* use.
- **The acetylation arm is not a proteasome strategy and is not treated as one here.** It is a
  different mechanism on the same protein — post-translational stability rather than bulk
  turnover — and it is developed as an **independent axis** in §5.3. It is separated out precisely
  because it does not inherit bortezomib's toxicity ceiling, which is the whole reason B sits
  last.

### 5.3 Hypothesis C — BubR1 acetylation-stability: an independent axis that multiplies A

**Why this is a separate axis and not a footnote to B.** Hypotheses A and B act on *different
objects*: translation termination (A) and bulk proteostasis (B). This axis acts on **the fraction of the protein that survives**. That makes
it the only lever in the report that **composes multiplicatively with readthrough**: A sets
how much BubR1 is *synthesised*, this sets how much is *retained*, and the two multiply. It is also
the only arm with a nutraceutical-grade safety ceiling for a child — no cytotoxicity, no
genotoxicity, no DSB.

**The biology.** BubR1 abundance is set post-translationally by acetylation. **SIRT2**, an
NAD⁺-dependent deacetylase, acts on BubR1 at two independently mapped sites:

- **K250** — deacetylation controls APC/C-dependent proteolysis and the timing of anaphase entry
  (PMID 25285631). This is the **human**-mapped site.
- **K668** — SIRT2 maintains the deacetylated state that the acetyltransferase **CBP** reverses;
  deacetylation *increases* BubR1 abundance. In that study, **NMN treatment raised BubR1 levels in
  vivo**, and SIRT2 overexpression extended median lifespan specifically in BubR1-hypomorphic mice
  (PMID 24825348). This is the **mouse** work.

Two results that read as objections in fact sharpen the case. SIRT2 overexpression had **no**
effect on lifespan in non-hypomorphic wild-type mice (PMID 38009412) — so the benefit is
**specific to animals under-expressing BubR1**, which is precisely this genotype, not a general
anti-ageing claim. And the direction of the two enzymes means the axis is falsifiable: SIRT2
activation and CBP inhibition should push the *same* way (raise BubR1), so an experiment in which
only one of them moves BubR1 abundance refutes the model rather than half-confirming it.

**Two ways to push the axis.**

- **Raise the substrate — NAD⁺ precursors: NMN, nicotinamide riboside, niacin.** The only arm in
  this report with decades of human exposure. **Caveat:** niacin is an approved drug; **NMN and NR
  are supplements, not approved medicines**, and their NAD⁺-raising effect is modest and
  variable. No human study has measured BubR1 abundance under any of them.
- **Remove the brake — p300/CBP inhibition: A-485 and CCS1477.** **A-485** is the selective
  *catalytic* p300/CBP inhibitor (PMID 28953875) but is a tool compound never dosed in patients;
  **CCS1477** is the *bromodomain* inhibitor actually in Phase 1/2 (prostate cancer). **Caveat:**
  bromodomain inhibition is not the same pharmacological event as catalytic inhibition, and no
  study has measured BubR1 acetylation under either agent.
  The mechanism is inferred, not shown.

**The allele-specificity warning, and it is the reason this axis cannot be prescribed blind.**
**K668 (668) sits N-terminal to the PTC (737)**, so it is present in the truncated allele-A
fragment. A stability-raising intervention applied without allele resolution could therefore
**stabilise the 737-aa fragment** — the dominant-negative risk that makes SMG1 inhibition alone
unacceptable (§5.1). The axis must be developed **against allele B**, with allele-A fragment
abundance measured as a safety readout, not assumed allele-agnostic.

**Two further gaps, stated plainly.** First, the species problem: human and mouse BubR1 share
almost no identity across residues 640–700 (3/61 identical), mouse residue 668 is **alanine**, and
mouse BubR1 is 1052 aa versus human 1050 — so the mouse K668 numbering does **not** transfer to
human K668 automatically, and the human acetylation map beyond K250 is incomplete. Second, none of
this has been tested in MVA cells: whether the missense allele is even subject to the same
acetylation-dependent turnover as wild-type BubR1 is unmeasured.

**The experiment that decides it, and it doubles as the A test.** A **2 × 2 design** in
patient-derived fibroblasts or iPSCs: ± readthrough agent × ± NAD⁺ precursor, with readouts
(a) full-length BubR1 protein, (b) **allele-A fragment abundance** (safety), (c) aneuploidy rate
after divisions, and (d) BubR1 acetylation by immunoprecipitation. If both arms raise BubR1 and the
combination is multiplicative, this is the cheapest two-hypothesis win available. If the fragment rises
alongside full-length protein, the axis is contraindicated until readthrough is established.

### 5.4 Investigational delivery mode: the aneuploid endolysosomal phenotype as a targeting handle

**A delivery mode, not a drug.** Nothing below is therapeutic on its own and none of it is
tested — it is a way to make A–B selective, and it stands or falls on one measurable ratio.

**The problem this solves.** A–B share one weakness: none of the agents is selective for the
child's diseased cells. Readthrough drugs act on every ribosome and every PTC in the body;
proteasome inhibitors are cytotoxic by design. Since there is no BubR1-specific
ligand to attach to (§4), the only remaining source of selectivity is **a property of the
aneuploid cell itself**.

**That property is the endolysosomal system — and, critically, the route into it changes.**
Aneuploid mammalian cells accumulate autophagosomal cargo — protein aggregates — inside lysosomes
that fail to degrade it, which triggers a **TFEB-mediated lysosomal stress response**: TFEB moves
into the nucleus and transcriptionally expands the autophagy/lysosomal gene programme
(PMID 26404941, 26571033). Constitutive trisomies show the same signature — upregulation of
autophagy and lysosomal degradation pathways, with constitutive TFEB nuclear localisation
(PMID 34215848) — and aneuploidy-driven rewiring of endo-lysosomal proteostasis is recognised as a
general feature (PMID 27896217).

**There is also direct evidence that the endocytic *route* shifts, not only the compartment.**
In human aneuploid blastocysts, RNA sequencing shows thousands of differentially expressed genes
clustered in **endocytosis pathways** (PSD3, ARAP2, SMAP2, CBLC, AGAP1, SH3GLB1), and fluorescent
probe measurements show the direction of the shift: **clathrin-dependent endocytosis is
reduced and clathrin-independent endocytosis is increased** in aneuploid versus diploid embryos
(PMID 41870862, P < 0.05). The authors attribute the shift to altered intracellular osmotic
pressure from an imbalanced amino-acid pool. **This is the observation that makes a delivery
strategy plausible at all**: an aneuploid cell is not merely a cell with more lysosomes, it is a
cell whose plasma membrane preferentially uses a *different* uptake route.

**The proposal: make the drug depend on the *clathrin-independent* route.** Do not merely make
the agent endocytosis-dependent — make it depend on the route that the aneuploid cell prefers and
the euploid cell does not. Concretely, carriers restricted to clathrin-independent uptake:
**macropinocytic** cargo (albumin conjugates, arginine-rich cell-penetrating peptides, particles
above the ~200 nm clathrin ceiling), **caveolae/lipid-raft**-targeted liposomes and raft-anchored
ligands, and **CLIC/GEEC** substrates — with an acid-labile or cathepsin-cleavable linker so the
payload is released only after it reaches the endolysosome. A clathrin-dependent carrier (the
transferrin-receptor route, most antibody–drug conjugates, most small anionic nanoparticles) would
be the *wrong* choice here: it targets the pathway aneuploid cells are losing. Three things
follow:

1. **Higher intracellular dose where the target is.** If aneuploid cells take up and retain more
   endocytic cargo and have greater lysosomal mass, an endocytosis-gated drug concentrates in
   exactly the cells that carry the extra chromosomes — readthrough in the cells that need
   BubR1.
2. **A lower systemic dose.** Which matters specifically here: the toxicity ceiling of the
   aminoglycoside/ELX-02 series is what limits them in a child (§5.1).
   Targeting converts a systemic-dose problem into a delivery problem.
3. **Selectivity without a BubR1 ligand** — the thing §4 says we cannot have.

**What is not yet established, stated plainly.** The route shift is measured in **human
embryos**, not in the child's somatic tissue, and it is a **single study**: whether it holds in
patient fibroblasts or iPSCs derived from this genotype is untested. Total endocytic *capacity*
is not quantified either — "more clathrin-independent" does not mean "more cargo overall".
Three further failure modes must be excluded before it is worth formulating anything: net
delivery may be **lower** despite the route shift; aneuploid cells may have **impaired**
degradative release of cargo, so a cleavable linker may never fire; lysosomal membrane
permeabilisation may be a toxicity route rather than a targeting route; and a constitutively
active TFEB/autophagy programme may make these cells *more* resistant to lysosomal death, not
less.

**The experiment that decides it is cheap and uses cells the programme needs anyway.** In
patient-derived fibroblasts or iPSCs versus euploid controls — and in aneuploid versus euploid
subclones sorted from the patient's own culture — measure, **separating the routes**:

- **clathrin-dependent** uptake: transferrin-AF647; blocked by chlorpromazine or hypertonic sucrose;
- **clathrin-independent** uptake: dextran-FITC and lactosylceramide-BSA (CLIC/GEEC), albumin or
  high-molecular-weight dextran ± EIPA (macropinocytosis), cholera-toxin B (raft/caveolae);
- plus **lysosomal mass and acidity** (LysoTracker, LAMP1, TFEB localisation) and net retention.

If the clathrin-independent : clathrin-dependent ratio is flat, this section is dead and A–B
stand as systemic-dose proposals. If the ratio shifts the way the embryo data predict,
delivery is the cheapest available selectivity in the whole programme — and it is the same cell
panel that measures residual BubR1 (§5.6).

### 5.5 Investigational therapeutic: base and prime editing — the definitive future fix

The three pharmacological hypotheses all *work around* the genotype: readthrough makes the ribosome misread the
stop; proteasome inhibition manages the consequences. **Only editing changes the
allele itself.** Because the mosaicism in MVA is of *aneuploidy*, not of genotype, **every cell
carries both variants** — there is no subpopulation to find, and under the candidate
compound-heterozygous model no wild-type allele to protect.
Any corrected cell gains function and there is no dominant-negative to preserve, which makes the
genetics unusually favourable. The barriers are delivery, timing and safety, not target
selection.

| Target | Codon change | Reversion needed | Editor that can do it |
|---|---|---|---|
| **Allele A** (p.Leu737Ter) | `TTA` (Leu) → `TGA` (stop) | **TGA→TGG (Trp)** — a single A→G — or exact reversion to `TTA` | **Adenine base editor (ABE)** for TGA→TGG; **prime editing** for exact Leu reversion |
| **Allele B** (p.Asn1002Lys) | `AAT` (Asn) → `AAG` (Lys) | **G→T** (a transversion) | **Prime editing only** — no base editor can convert G→T |

**Allele A is the approachable one.** An ABE converts the opal codon to tryptophan with **no
double-strand break and no HDR** — exactly what a checkpoint-defective genome needs. Residue
737 lies N-terminal to the pseudokinase domain (766–1050) and Leu→Trp is a conservative
hydrophobic substitution, so the corrected protein should retain scaffolding function — but
this must be **confirmed, not assumed**. ABE8e-class editors have already corrected a human
nonsense mutation precisely and efficiently (PMID 36385635), and base editors have now been
delivered non-virally to correct nonsense variants in vivo (PMID 42276271) and to the CNS
(PMID 42026681) — the two barriers that matter here.

**Allele B needs prime editing** and is the harder prize: prime editing installs arbitrary
substitutions without DSBs, but efficiency in post-mitotic tissue is the open question. A third
route worth naming is **RNA-level editing (ADAR-based)**, which leaves the genome untouched and
is re-dosable — attractive where permanent editing carries permanent risk. It is listed as an
option, not yet developed for this target.

Figure 7 (above, §5.1) is the codon-level view: the same two codons that set the readthrough
problem set the editing problem.

**Why it is still not a first-line recommendation today:** (i) **body-wide delivery** to a child's
brain, gut and marrow has no clinical route; (ii) the **developmental window** in which most of
the phenotype is laid down largely precedes diagnosis; (iii) **per-cell efficiency** must clear
the same ~50% BubR1 threshold, tissue by tissue; (iv) it is **not drug repurposing**, and
carries a separate regulatory and ethical path.

**The experiment that would justify it:** take patient-derived iPSCs, apply ABE (TGA→TGG) and
prime editing (allele B) in parallel, and measure (a) % full-length BubR1 restored and (b)
aneuploidy rate under spindle stress. That is a bench-scale, falsifiable test of the definitive
fix — and it uses the same readout A needs as its baseline.

### 5.6 The single experiment that decides everything

**Measure BubR1 in patient-derived cells, as % of control.** That number sets the readthrough
target (`50% − measured%`), tells us whether the programme is feasible, and serves as the
baseline for the readthrough screen and the editing benchmark. Until it is measured, this is a
well-reasoned, falsifiable hypothesis — not a treatment.

Run it as a **2 × 2**: ± readthrough agent × ± NAD⁺ precursor, with full-length BubR1,
**allele-A fragment abundance**, and aneuploidy rate after divisions as readouts (§5.3). One
experiment then prices two hypotheses and tests whether they compose, which is the only way to know
whether the combination is worth pursuing before any of it reaches an animal.

---

## 6. Diagnostics

**This section proposes tests, not treatments.** A molecular diagnosis is what would let a case
like this one be recognised before it is a crisis, and a biomarker of the aneuploid-cell burden is
what would let any of §5 be monitored. Two complementary tools are proposed: a **gene panel** that
asks *which gene* (§6.1), and **SASP markers in cord blood** that ask *how much aneuploid-cell
stress is already present at birth* (§6.2). Both are research-use only and neither has been
validated.

### 6.1 The ~120-gene MVA panel as a diagnostic tool

We curated a **123-gene MVA panel** (the "~120-gene panel") — 10 with direct MVA evidence (*BUB1B*, *CEP57*, *TRIP13*,
*MAD1L1*, *SMC5*, *TRIM37* and others), 95 curated mitotic-machinery genes (SAC, kinetochore,
centrosome, cohesin, SMC5/6), and 18 CIN-adjacent context genes — from GenCC, ClinGen, PanelApp,
OMIM, Orphadata, G2P, HGNC and Complex Portal, with a GRCh38 BED and dated source versions
(`analysis/2026-10-23_mva-gene-panel-sources/`). Run phenotype-agnostically over this child's WGS it
recovered the *BUB1B* pair as the only genotype-plausible cause and found no second biallelic gene with
direct MVA evidence — an independent negative, not a re-ranking
(`analysis/2026-10-23_panel-vs-exomiser-overlay/`).

Its natural use is **a new index case, prenatally**. When a fetus already shows the MVA phenotype
— growth restriction, microcephaly, mosaic aneuploidy on karyotype or interphase FISH — amniocyte
or CVS DNA can be run through the same 123 genes, and unlike a singleton WGS it answers the
question the family actually asks: *which gene, and therefore what recurrence risk*. Two things it
cannot do: it cannot establish the MVA phenotype itself, which stays cytogenetic; and it cannot make
a VUS diagnostic — allele B in this child would still come back as one. Once a family's alleles are
known, a targeted assay for those two variants is faster, cheaper and more sensitive than a panel,
and is the right test for a subsequent pregnancy.

Three caveats travel with any such claim. Performance on prenatal samples is **untested** — low
input, whole-genome amplification and maternal cell contamination all degrade it, and contamination
must be excluded explicitly. The unsolved fraction of MVA means a **negative panel excludes
nothing**. And clinical use requires a validated assay and accreditation, which this is not: it is
diagnostic *content*, research-use only.

### 6.2 Aneuploid cells, SASP, and markers measurable in cord blood

**The rationale.** Aneuploid cells frequently arrest and adopt a **senescent** phenotype, and
senescent cells secrete the **senescence-associated secretory phenotype (SASP)** — a mix of
cytokines, chemokines and proteases that drives chronic inflammation. Aneuploidy-driven senescence
is documented: trisomy 21 induces senescence that disrupts the nuclear architecture of neural
progenitors (PMID 34995493), and clearance of senescent cells has been proposed for aneuploidy-prone
tissue (PMID 30795536). Constitutive trisomies also show an innate-immune / cGAS-STING signature
(PMID 34215848). In MVA the same logic predicts a **secreted, measurable footprint of the aneuploid
cell burden** — which is what makes it a candidate *diagnostic* readout rather than a therapeutic
target. No MVA or *BUB1B* data exist for it; this is extrapolation from trisomy-21 and ageing biology.

**Why cord blood.** It is the one fluid that can be sampled at birth without invasive procedures,
before treatment and before secondary illness confounds the picture, and it allows a newborn with
prenatal growth restriction or microcephaly to be flagged early. It is the postnatal counterpart to
the prenatal use of the panel in §6.1.

**What to measure (a proposed panel, not a validated one).**

- **Core SASP cytokines/chemokines:** IL-6, IL-8 (CXCL8), IL-1β, MCP-1 (CCL2), GDF-15.
- **SASP matrix/remodelling factors:** MMP-3, PAI-1 (SERPINE1), ICAM-1.
- **Innate-immune read-out:** interferon-stimulated signature, as a marker of the cGAS-STING
  response to chromosome mis-segregation.
- **A direct aneuploidy anchor, not a surrogate:** interphase FISH or karyotype on cord-blood
  lymphocytes, so that any secreted marker is read against the mosaic aneuploid fraction actually
  present.
- **Genotype in parallel:** the panel (§6.1) run on the same sample, so a raised marker profile is
  interpreted against a known BUB1B/MVA genotype, not in isolation.

**What this cannot do, stated plainly.**

- **Specificity is poor.** IL-6, IL-8 and CRP-type markers rise sharply with chorioamnionitis,
  prematurity, birth asphyxia and neonatal sepsis; a raised cord-blood SASP profile is not
  diagnostic of aneuploidy, and a normal one does not exclude it. It can support, never establish.
- **Cord blood may under-report.** The aneuploid cells may sit in tissues (brain, gut, marrow
  niches) that do not shed into the circulation, and the senescent fraction in a neonate may be
  small.
- **No reference ranges or cut-offs exist** for MVA; a case–control cohort of affected and euploid
  neonates would be required before any threshold could be set. Mosaicism is rare, so that cohort
  is hard to assemble.
- **It is a monitoring and triage tool.** If SASP markers do track aneuploid burden, the same
  panel could serve as a pharmacodynamic readout for §5 — but that is a second hypothesis resting
  on the first.

A preliminary maternal-serum proteomics report on aneuploid pregnancies (listed under *Sources
cited without a PMID* in §8; n = 5, hypothesis-generating only) points the same way — that
aneuploidy leaves a protein footprint in a circulating fluid — but is not evidence for the cord-blood
panel proposed here.

**The experiment that decides it** reuses cells the programme needs anyway: measure the §6.2 marker
set in conditioned medium from patient-derived fibroblasts or iPSCs versus euploid controls, and in
aneuploid versus euploid subclones sorted from the patient's own culture. If the secreted profile
does not separate them *in vitro*, there is no reason to take it to cord blood.

---

## 7. Acknowledgements

**Team SherlockMVA**

- **Priyanka Rajan**
- **Soumya Tanwar**
- **Shreetoma Chakraborty**

The team collectively carried out proband ascertainment and variant interpretation,
structural bioinformatics (model preparation, docking, fpocket), molecular dynamics
(setup, production runs, analysis), literature curation with PMID-level verification of
every external claim, and the preparation of this report and its figures.

We thank the family for entrusting us with the proband's data, and the clinical genetics and
care teams who supported the diagnosis. Computational work was performed on shared local
infrastructure; all analyses are archived and re-runnable
(`analysis/2026-10-20_bub1b-missense-functional-prediction/`,
`analysis/2026-10-23_bub1b-virtual-docking/`, `analysis/2026-10-23_bub1b-md-stability/`,
`analysis/2026-10-23_track2-lit-scan/`, `analysis/2026-10-24_track2-reconciliation/`).

## 8. References

Numbered in order of first citation in the text; the PMID is repeated inline at each point of
use, so either the number here or the PMID in the sentence leads to the same source. All
citations were resolved against PubMed E-utilities on 2026-10-24 (see
`analysis/2026-10-24_track2-reconciliation/results/verification.txt`).

1. Gama Braga L, Cisneros AF, Mathieu MM, et al. BUBR1 Pseudokinase Domain Promotes Kinetochore
   PP2A-B56 Recruitment, Spindle Checkpoint Silencing, and Chromosome Alignment. *Cell Rep*.
   2020;33(7):108397. doi:10.1016/j.celrep.2020.108397. PMID: 33207204.
2. Crawford DK, Alroy I, Sharpe N, et al. ELX-02 Generates Protein via Premature Stop Codon
   Read-Through without Inducing Native Stop Codon Read-Through Proteins. *J Pharmacol Exp Ther*.
   2020;374(2):264-272. doi:10.1124/jpet.120.265595. PMID: 32376628.
3. Destoop M, Ramalho AS, Silva IAL, et al. A functional comparison of readthrough agent ELX-02
   across a wide range of nonsense CFTR variants. *J Cyst Fibros*. 2026. pii:
   S1569-1993(26)03738-0. doi:10.1016/j.jcf.2026.09.004. PMID: 42728176.
4. Berger J, Li M, Berger S, et al. Effect of Ataluren on dystrophin mutations. *J Cell Mol Med*.
   2020;24(12):6680-6689. doi:10.1111/jcmm.15319. PMID: 32343037.
5. Mosallaei D, Hao M, Antaya RJ, et al. Molecular and Clinical Outcomes After Intravenous
   Gentamicin Treatment for Patients With Junctional Epidermolysis Bullosa Caused by Nonsense
   Variants. *JAMA Dermatol*. 2022;158(4):366-374. doi:10.1001/jamadermatol.2021.5992.
   PMID: 35234826.
6. Lebeda D, Smoljanow D, Hofhuis J, et al. Release factor 3 degradation enhances
   aminoglycoside-mediated translational readthrough at MECP2 nonsense mutations. *J Mol Med
   (Berl)*. 2026;104(1):123. doi:10.1007/s00109-026-02728-z. PMID: 42907672.
7. Chen J, Zhou M, Dong S, et al. Nonviral delivery of chemically modified tRNA rescues nonsense
   mutations in cystic fibrosis. *Science*. 2026;393(6814):eaeb0054. doi:10.1126/science.aeb0054.
   PMID: 42658930.
8. Bohers E, Sarafan-Vasseur N, Drouet A, et al. Gradual reduction of BUBR1 protein levels
   results in premature sister-chromatid separation then in aneuploidy. *Hum Genet*.
   2008;124(5):473-478. doi:10.1007/s00439-008-0572-y. PMID: 18932004.
9. Licciardi F, Lhakhang T, Kramer YG, et al. Human blastocysts of normal and abnormal
    karyotypes display distinct transcriptome profiles. *Sci Rep*. 2018;8(1):14906.
    doi:10.1038/s41598-018-33279-0. PMID: 30297919.
10. Tsai HJ, Nelliat AR, Choudhury MI, et al. Hypo-osmotic-like stress underlies general cellular
    defects of aneuploidy. *Nature*. 2019;570(7759):117-121. doi:10.1038/s41586-019-1187-2.
    PMID: 31068692.
11. Ippolito MR, Zerbib J, Eliezer Y, et al. Increased RNA and Protein Degradation Is Required
    for Counteracting Transcriptional Burden and Proteotoxic Stress in Human Aneuploid Cells.
    *Cancer Discov*. 2024;14(12):2532-2553. doi:10.1158/2159-8290.CD-23-0309. PMID: 39247952.
12. Suematsu T, Li Y, Kojima H, et al. Deacetylation of the mitotic checkpoint protein BubR1 at
    lysine 250 by SIRT2 and subsequent effects on BubR1 degradation during the
    prometaphase/anaphase transition. *Biochem Biophys Res Commun*. 2014;453(3):588-594.
    doi:10.1016/j.bbrc.2014.09.128. PMID: 25285631.
13. North BJ, Rosenberg MA, Jeganathan KB, et al. SIRT2 induces the checkpoint kinase BubR1 to
    increase lifespan. *EMBO J*. 2014;33(13):1438-1453. doi:10.15252/embj.201386907.
    PMID: 24825348.
14. Wu LE, Fiveash CE, Bentley NL, et al. SIRT2 transgenic over-expression does not impact
    lifespan in mice. *Aging Cell*. 2023;22(12):e14027. doi:10.1111/acel.14027. PMID: 38009412.
15. Lasko LM, Jakob CG, Edalji RP, et al. Discovery of a selective catalytic p300/CBP inhibitor
    that targets lineage-specific tumours. *Nature*. 2017;550(7674):128-132.
    doi:10.1038/nature24028. PMID: 28953875.
16. Santaguida S, Vasile E, White E, et al. Aneuploidy-induced cellular stresses limit autophagic
    degradation. *Genes Dev*. 2015;29(19):2010-2021. doi:10.1101/gad.269118.115. PMID: 26404941.
17. Krivega M, Stiefel CM, Karbassi S, et al. Genotoxic stress in constitutive trisomies induces
    autophagy and the innate immune response via the cGAS-STING pathway. *Commun Biol*.
    2021;4(1):831. doi:10.1038/s42003-021-02278-9. PMID: 34215848.
18. Demirsoy S, Martin S, Maes H, et al. Adapt, Recycle, and Move on: Proteostasis and
    Trafficking Mechanisms in Melanoma. *Front Oncol*. 2016;6:240. doi:10.3389/fonc.2016.00240.
    PMID: 27896217.
19. Duan R, Liang X, Yang Y, et al. Abnormal Endocytosis Resulting from Imbalanced Proteomic
    Doses in Human Aneuploid Embryos. *Reprod Sci*. 2026;33(4):774-785.
    doi:10.1007/s43032-026-02074-y. PMID: 41870862.
20. Sheriff A, Guri I, Zebrowska P, et al. ABE8e adenine base editor precisely and efficiently
    corrects a recurrent COL7A1 nonsense mutation. *Sci Rep*. 2022;12(1):19643.
    doi:10.1038/s41598-022-24184-8. PMID: 36385635.
21. Hiya-Kawaguchi U, Kashiwakura Y, Baatartsogt N, et al. Nonviral delivery of a base editor
    enables personalized correction of hemophilia B nonsense variants in a mouse model. *J
    Thromb Haemost*. 2026;24(9):3159-3171. doi:10.1016/j.jtha.2026.06.003. PMID: 42276271.
22. Nam BG, Seo JH, Hong SA, et al. In vivo adenine base editing of mutant Galc gene ameliorates
    Krabbe disease progression. *Genome Med*. 2026;18(1):46. doi:10.1186/s13073-026-01620-2.
    PMID: 42026681.
23. Meharena HS, Marco A, Dileep V, et al. Down-syndrome-induced senescence disrupts the nuclear
    architecture of neural progenitors. *Cell Stem Cell*. 2022;29(1):116-130.e7.
    doi:10.1016/j.stem.2021.12.002. PMID: 34995493.
24. Melo Pereira S, Ribeiro R, Logarinho E. Approaches towards Longevity: Reprogramming,
    Senolysis, and Improved Mitotic Competence as Anti-Aging Therapies. *Int J Mol Sci*.
    2019;20(4):938. doi:10.3390/ijms20040938. PMID: 30795536.

**Sources cited in the text without a PMID** (grey literature / databases, listed for
completeness; not peer-reviewed): OMIM 257300 (MVA1); gnomAD v4.1 allele frequency; AlphaMissense,
ESM-1v and REVEL in-silico predictions (interim evidence, `evidence/`); and the preliminary
maternal-serum proteomics report on aneuploid pregnancies (López-Uriarte GA et al., *Medical
Research Archives* 2026;14(1), doi:10.18103/mra.v14i1.7199 — not MEDLINE-indexed, n = 5,
hypothesis-generating only; see `analysis/2026-10-24_track2-reconciliation/results/diagnostic_paper.md`).

---

*This report is a hypothesis for further investigation, not evidence that any medicine works.
It is not medical advice and must not be used to justify off-label treatment of this child
outside a supervised clinical trial or compassionate-use protocol.*
