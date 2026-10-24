// =====================================================================
// BUB1B N1002K — unified reproducible pipeline (Nextflow DSL2)
//
// Stage 1  structural substrate / druggability   (env: bubdock)
//   FETCH_AF -> EXTRACT_KINASE -> MUTATE
//            -> FPOCKET(full, WT-kinase, mut-kinase) -> SUMMARIZE
//            -> POCKET_SURVEY (97-pocket annotated survey)
//
// Stage 2  WT vs N1002K molecular dynamics         (env: bubmd)
//   PREP_SYSTEM(WT, mut) -> MD(WT, mut, HMR/4fs)
//                        -> ANALYZE (pseudokinase + stability reports)
//
// Rationale and interpretation: see README.md and the parent analysis dirs.
// Runs with the `local` executor (no Docker/Singularity on this host) and
// micromamba envs. Launch via:  bash run_pipeline.sh [ -resume ]
// =====================================================================
nextflow.enable.dsl = 2

// ---------------------------------------------------------------------
// Stage 1: structure & druggability
// ---------------------------------------------------------------------
process FETCH_AF {
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/01_structure", mode: 'copy'
    input:
    val acc
    output:
    path "${acc}.pdb",       emit: pdb
    path "${acc}_plddt.tsv", emit: plddt
    stub:
    """
    touch ${acc}.pdb ${acc}_plddt.tsv
    """
    script:
    """
    micromamba run -n ${params.env_struct} python ${projectDir}/bin/fetch_af.py ${acc} ${acc}.pdb ${acc}_plddt.tsv
    """
}

process EXTRACT_KINASE {
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/02_regions", mode: 'copy'
    input:
    path pdb
    output:
    path "bub1b_kinase.pdb", emit: pdb
    stub:
    """
    touch bub1b_kinase.pdb
    """
    script:
    """
    micromamba run -n ${params.env_struct} python ${projectDir}/bin/extract_region.py ${pdb} bub1b_kinase.pdb ${params.kstart} ${params.kend}
    """
}

process MUTATE {
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/02_regions", mode: 'copy'
    input:
    path pdb
    output:
    path "bub1b_kinase_N1002K.pdb", emit: pdb
    stub:
    """
    touch bub1b_kinase_N1002K.pdb
    """
    script:
    """
    micromamba run -n ${params.env_struct} pymol -cq ${projectDir}/bin/mutate.pml -- ${pdb} bub1b_kinase_N1002K.pdb ${params.mut_resid} ${params.mut_aa}
    """
}

process FPOCKET {
    tag "$region"
    cpus 2
    memory '3 GB'
    publishDir "${params.outdir}/03_druggability", mode: 'copy'
    input:
    tuple val(region), path(pdb)
    output:
    path "pocket_summary_${region}.tsv", emit: tsv
    stub:
    """
    touch pocket_summary_${region}.tsv
    """
    script:
    """
    micromamba run -n ${params.env_struct} fpocket -f ${pdb}
    micromamba run -n ${params.env_struct} python ${projectDir}/bin/fpocket_summary.py ${pdb} ${region} pocket_summary_${region}.tsv ${params.mut_resid}
    """
}

process POCKET_SURVEY {
    cpus 2
    memory '3 GB'
    publishDir "${params.outdir}/04_pocket_survey", mode: 'copy'
    input:
    path pdb
    path plddt
    output:
    path "pockets_annotated.tsv", emit: tsv
    stub:
    """
    touch pockets_annotated.tsv
    """
    script:
    """
    micromamba run -n ${params.env_struct} python ${projectDir}/bin/pocket_survey.py ${pdb} ${plddt} .
    """
}

process SUMMARIZE {
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/03_druggability", mode: 'copy'
    input:
    path tsvs
    output:
    path "fpocket_all_regions.tsv", emit: all
    stub:
    """
    touch fpocket_all_regions.tsv
    """
    script:
    """
    micromamba run -n ${params.env_struct} python ${projectDir}/bin/combine.py fpocket_all_regions.tsv ${tsvs}
    """
}

// ---------------------------------------------------------------------
// Stage 2: WT vs N1002K molecular dynamics
// ---------------------------------------------------------------------
process PREP_SYSTEM {
    tag "$label"
    cpus 1
    memory '3 GB'
    publishDir "${params.outdir}/05_md_prep", mode: 'copy'
    input:
    tuple val(label), path(pdb)
    output:
    tuple val(label), path("kinase_${label}.pdb"), emit: ready
    stub:
    """
    touch kinase_${label}.pdb
    """
    script:
    """
    micromamba run -n ${params.env_md} python ${projectDir}/bin/prep_system.py --inp ${pdb} --out kinase_${label}.pdb --ph 7.0
    """
}

process MD {
    tag "$label"
    cpus params.md_threads
    memory '6 GB'
    publishDir "${params.outdir}/06_md", mode: 'copy'
    input:
    tuple val(label), path(pdb)
    output:
    tuple val(label), path("kinase_${label}.dcd"), path("kinase_${label}.final.pdb"), path("kinase_${label}.log.csv"), emit: traj
    stub:
    """
    touch kinase_${label}.dcd kinase_${label}.final.pdb kinase_${label}.log.csv
    """
    script:
    """
    micromamba run -n ${params.env_md} python ${projectDir}/bin/md_run.py \
        --pdb ${pdb} --outprefix kinase_${label} \
        --ns ${params.md_ns} --equil-ns ${params.md_equil_ns} \
        --threads ${params.md_threads} --save-ps ${params.save_ps} --hmr
    """
}

// ANALYSIS reads the published trajectories from ${outdir}/06_md (absolute),
// and is gated on all MD tasks finishing via the collected `traj` channel.
process ANALYZE {
    cpus 2
    memory '4 GB'
    publishDir "${params.outdir}/07_md_analysis", mode: 'copy'
    input:
    val done
    output:
    path '*', emit: analysis
    stub:
    """
    touch pseudokinase_summary.tsv pseudokinase_verdict.txt summary_stability.tsv verdict.txt
    """
    script:
    def md_dir = "${projectDir}/${params.outdir}/06_md"
    """
    micromamba run -n ${params.env_md} python ${projectDir}/bin/analyze_pseudokinase.py --dir ${md_dir} --out . --target ${params.mut_resid}
    micromamba run -n ${params.env_md} python ${projectDir}/bin/analyze_stability.py --dir ${md_dir} --out . --target ${params.mut_resid}
    """
}

// ---------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------
workflow {
    FETCH_AF(params.acc)

    EXTRACT_KINASE(FETCH_AF.out.pdb)
    MUTATE(EXTRACT_KINASE.out.pdb)

    // --- Stage 1: druggability survey ---
    ch_full   = FETCH_AF.out.pdb.map        { p -> tuple("bub1b_full", p) }
    ch_kinwt  = EXTRACT_KINASE.out.pdb.map  { p -> tuple("bub1b_kinase_WT", p) }
    ch_kinmut = MUTATE.out.pdb.map          { p -> tuple("bub1b_kinase_N1002K", p) }

    FPOCKET(ch_full.mix(ch_kinwt).mix(ch_kinmut))
    SUMMARIZE(FPOCKET.out.tsv.collect())
    POCKET_SURVEY(FETCH_AF.out.pdb, FETCH_AF.out.plddt)

    // --- Stage 2: MD stability (optional) ---
    if (params.run_md) {
        ch_wt  = EXTRACT_KINASE.out.pdb.map { p -> tuple("WT", p) }
        ch_mut = MUTATE.out.pdb.map         { p -> tuple("N1002K", p) }

        PREP_SYSTEM(ch_wt.mix(ch_mut))
        MD(PREP_SYSTEM.out.ready)

        // gate ANALYZE on completion of all MD trajectories
        ANALYZE(MD.out.traj.collect())
    }
}
