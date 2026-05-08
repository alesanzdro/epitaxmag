#!/usr/bin/env nextflow

/*
 * ══════════════════════════════════════════════════════════════════
 *  EpiTaxMAG — Nanopore metagenomics pipeline
 *  TAX phase  (EpiTax):    QC · trim · filter · taxonomic profiling
 *  MAG phase  (EpiTaxMAG): assembly · binning · annotation · AMR
 *  EPIMOL — FISABIO
 * ══════════════════════════════════════════════════════════════════
 *
 *  USAGE:
 *    # Taxonomic profiling only (EpiTax)
 *    nextflow run main.nf --input data/260225_EPIM230 -profile gru
 *
 *    # Full pipeline: profiling + MAG recovery (EpiTaxMAG)
 *    nextflow run main.nf --input data/260225_EPIM230 -profile gru --run_assembly
 *
 *    # MAG only (requires previously filtered reads)
 *    nextflow run main.nf --input data/260225_EPIM230 -profile gru --run_tax false --run_assembly
 *
 *    # Resume after interruption
 *    nextflow run main.nf --input data/260225_EPIM230 -profile gru -resume
 */

nextflow.enable.dsl=2


// ══════════════════════════════════════════════════════════════════
// IMPORTS — subworkflows
// ══════════════════════════════════════════════════════════════════

include { SETUP_DB as EPITAXMAG_SETUP } from './workflow/setup_db'
include { TAX as EPITAXMAG_TAX } from './workflow/tax/main'
include { MAG as EPITAXMAG_MAG } from './workflow/mag/main'


// ══════════════════════════════════════════════════════════════════
// HELP
// ══════════════════════════════════════════════════════════════════

if (params.help) {
    log.info """
    EpiTaxMAG — Nanopore metagenomics pipeline
    EPIMOL — FISABIO
    ──────────────────────────────────────────

    USAGE:
      nextflow run main.nf --input <fastq_dir> -profile <profile> [options]

    MAIN OPTIONS:
      --input            Directory with .fastq.gz files (REQUIRED)
      --db_root          Database root [default: ${projectDir}/databases]
      --outdir           Output directory [default: results/<run_name>]
      --run_name         Run name [default: input directory name]

    PHASES:
      --run_tax          Run TAX phase [default: true]
      --run_assembly     Run MAG phase [default: false]

    READ FILTERING:
      --min_quality      Minimum Chopper quality [default: 10]
      --min_length       Minimum read length (bp) [default: 1000]
      --host_genome      Host reference genome (.fasta) for HOST_REMOVAL.
                         If omitted, host depletion is skipped.

    TAXONOMY:
      --kraken2_confidence  Kraken2 confidence [default: 0.1]
      --genome_size_mb      Average genome size in Mb [default: 4.5]
      --min_coverage_mag    Minimum theoretical coverage for MAGs [default: 30]

    PHENOTYPIC VERIFICATION (minimap2 vs type-strains):
      --phenotypic_targets  TSV with target organisms + NCBI accessions
                            (columns: organism, accession). Enabled by default.
                            Edit assets/phenotypic_targets/target_organisms.tsv
                            to add/remove organisms.
                            To disable: --phenotypic_targets false

    FASTQ SCREEN:
      --fastqscreen_conf    FastQ Screen config (.conf)

    DATABASES:
      --db_root             Database root [default: databases/]
      --skip_db_setup       Skip automatic DB download [default: false]

    MAG OPTIONS:
      --taxonomy_tool       MAG taxonomy classifier [default: skani]
                            Values:
                              skani    — fast ANI vs GTDB r226, ~16-24 GB RAM, ~36 GB DB
                              sourmash — k-mer sketches, ~5 GB RAM, ~3.7 GB DB
                              gtdbtk   — gold-standard, ≥120 GB RAM, ~60 GB DB
      --assemblies_dir      Reuse existing assemblies (skips MetaFlye)
      --filtered_dir        Reuse pre-filtered reads (skips TAX)
      --skip_gtdbtk         DEPRECATED — use --taxonomy_tool sourmash instead.
                            Kept as an alias that maps to 'sourmash'.

    PROFILES (recommended: -profile medium):
      -profile low          Modest workstation: ≥16 threads, ≥32 GB RAM
      -profile medium       Default tier:        ≥24 threads, ≥64 GB RAM
      -profile high         Server-class:        ≥64 threads, ≥128 GB RAM
      -profile singularity  Use Singularity containers (combinable: -profile medium,singularity)

      DEPRECATED ALIASES (still working): pc=low, censalud=medium, gru=high.

    EXAMPLES:
      # Taxonomic profiling with conda (default profile = medium)
      nextflow run main.nf --input data/EPIM230 -profile medium

      # Modest workstation (16 cores / 32 GB RAM)
      nextflow run main.nf --input data/EPIM230 -profile low

      # Server-class host with GTDB-Tk gold-standard taxonomy
      nextflow run main.nf --input data/EPIM230 -profile high --run_assembly --taxonomy_tool gtdbtk

      # With Singularity
      nextflow run main.nf --input data/EPIM230 -profile medium,singularity

      # With human host removal
      nextflow run main.nf --input data/EPIM230 -profile medium \\
          --host_genome /path/to/GRCh38.fna

      # Custom database root
      nextflow run main.nf --input data/EPIM230 -profile medium \\
          --db_root /disk/databases

      # Disable phenotypic verification
      nextflow run main.nf --input data/EPIM230 -profile medium \\
          --phenotypic_targets false

      # Resume after interruption
      nextflow run main.nf --input data/EPIM230 -profile medium -resume
    """.stripIndent()
    System.exit(0)
}


// ══════════════════════════════════════════════════════════════════
// VALIDATION
// ══════════════════════════════════════════════════════════════════

if (!params.input) {
    error "ERROR: --input is required. Usage: nextflow run main.nf --input /path/to/fastqs -profile <gru|pc|censalud>"
}

// ── Backwards-compatibility shim for the deprecated `--skip_gtdbtk` flag.
// If a user still passes it, route them onto Sourmash (which is what
// `--skip_gtdbtk` used to do) and warn loudly — the new selector is
// `--taxonomy_tool {skani|sourmash|gtdbtk}` (default: 'skani').
if (params.skip_gtdbtk && params.taxonomy_tool == 'skani') {
    log.warn "--skip_gtdbtk is deprecated; use --taxonomy_tool sourmash (or gtdbtk) instead. Falling back to 'sourmash'."
    params.taxonomy_tool = 'sourmash'
}

// Validate the selector early so misspellings fail fast.
def valid_tools = ['skani', 'sourmash', 'gtdbtk']
if (!(params.taxonomy_tool in valid_tools)) {
    error "Unknown --taxonomy_tool '${params.taxonomy_tool}'. Choose one of: ${valid_tools.join(', ')}."
}

// ── Profile validation and deprecation warnings ─────────────────────
// Users should pass `-profile {low|medium|high}`. The names `pc`,
// `censalud` and `gru` are deprecated aliases kept for backwards
// compatibility — they map to low/medium/high respectively.
def profile_aliases = [pc: 'low', censalud: 'medium', gru: 'high']
def active_profile  = (workflow.profile?.split(',') ?: [])
                      .findAll { it && it != 'standard' && it != 'singularity' }
                      .find { true }   // first matched

if (!active_profile) {
    log.warn """
    No size profile selected. Add `-profile medium` (recommended default) or
    one of: low (≥16t/32GB) | medium (≥24t/64GB) | high (≥64t/128GB).
    Falling back to whatever the global config provides — process resources
    may be undersized.
    """.stripIndent()
} else if (active_profile in profile_aliases) {
    log.warn "Profile '${active_profile}' is deprecated; use '${profile_aliases[active_profile]}' instead."
}

// GTDB-Tk requires ≥120 GB RAM. Reject early on profiles that cannot host it.
if (params.taxonomy_tool == 'gtdbtk' && active_profile && active_profile != 'high' && active_profile != 'gru') {
    error """
    --taxonomy_tool gtdbtk requires the 'high' profile (≥128 GB RAM, ≥64 CPUs)
    because pplacer alone needs ≥120 GB for the GTDB r232 reference data.
    Either:
      • Switch to the recommended classifier:  --taxonomy_tool skani  (default; ~24 GB RAM)
      • Use a lighter alternative:              --taxonomy_tool sourmash  (~8 GB RAM)
      • Run on a server-class host:             -profile high
    """.stripIndent()
}

def run_name    = params.run_name ?: file(params.input).name
def pipeline_id = params.run_assembly ? 'EpiTaxMAG' : 'EpiTax'

// ── TAX database checks (warn, do not abort — auto-setup will run) ─
def db_missing = []
if (params.run_tax) {
    def db = params.db_root
    def tax_checks = [
        ["${db}/k2_pluspf_20251015",                    'hash.k2d',  'Kraken2 PlusPF'],
        ["${db}/kaiju/kaiju_db_refseq_ref.fmi",         null,        'Kaiju FMI'],
        ["${db}/sylph/gtdb-r226-c200-dbv1.syldb",       null,        'Sylph GTDB'],
        ["${db}/kma_resfinder/resfinder_db/all.fsa",    null,        'ResFinder'],
    ]
    tax_checks.each { path, marker, name ->
        def check_path = marker ? "${path}/${marker}" : path
        if (!file(check_path).exists()) {
            db_missing << "  - ${name}: ${check_path}"
        }
    }
}
if (db_missing && params.skip_db_setup) {
    log.warn """
    ┌─────────────────────────────────────────────────────────────┐
    │  MISSING DATABASES                                          │
    └─────────────────────────────────────────────────────────────┘
    ${db_missing.join('\n    ')}

    --skip_db_setup is enabled, so no automatic download will be performed.
    Either remove --skip_db_setup, or download the missing files manually.
    """.stripIndent()
} else if (db_missing) {
    log.info """
    Missing databases will be auto-downloaded by SETUP_DB:
    ${db_missing.join('\n    ')}
    """.stripIndent()
}

log.info """
  =========================================
   E p i T a x M A G   v1.0
  =========================================
   Nanopore Metagenomics Pipeline
   FISABIO - EPIMOL
  -----------------------------------------

  Run:        ${run_name}
  Input:      ${params.input}
  Output:     ${params.outdir ?: "${projectDir}/results/${run_name}"}
  DB Root:    ${params.db_root}
  TAX phase:  ${params.run_tax ? 'YES' : 'NO'}
  MAG phase:  ${params.run_assembly ? 'YES' : 'NO'}
  Filters:    Q>=${params.min_quality}, len>=${params.min_length}bp
  MAG cov:    >=${params.min_coverage_mag}x (genome ${params.genome_size_mb} Mb)
""".stripIndent()


// ══════════════════════════════════════════════════════════════════
// INPUT CHANNEL — discover samples (skips files matching 'unclassified')
// ══════════════════════════════════════════════════════════════════

Channel
    .fromPath("${params.input}/*.fastq.gz")
    .filter { !it.name.contains('unclassified') }
    .map { file ->
        def sample = file.baseName.replaceAll(/\.fastq$/, '')
        tuple(sample, file)
    }
    .set { ch_raw_fastq }


// ══════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ══════════════════════════════════════════════════════════════════

workflow {

    // ── Setup: automatic database download ──────────────────────
    // SETUP_DB writes every required database into params.db_root via
    // `storeDir`, so it is idempotent across runs. We deliberately do
    // NOT chain SETUP_DB into TAX/MAG via a synthetic "ready" channel:
    // doing so would change the upstream hash signature on every run
    // and invalidate the resume cache for every TAX/MAG task.
    //
    // On a first run with empty databases there is a small race window
    // for SETUP_FASTQSCREEN (it has to download a 1.1 GB Zenodo panel
    // before TAX:FASTQSCREEN consumes the generated config). The TAX
    // process retries once on failure (errorStrategy retry); if both
    // attempts fall inside the download window, just relaunch with
    // `-resume` and Nextflow will pick up the cached setup outputs.
    if (!params.skip_db_setup) {
        EPITAXMAG_SETUP()
    }

    // ── TAX phase: taxonomic profiling (EpiTax) ─────────────────
    if (params.run_tax) {
        EPITAXMAG_TAX(ch_raw_fastq)
    }

    // ── MAG phase: assembly + MAG recovery (EpiTaxMAG) ──────────
    if (params.run_assembly) {
        if (params.run_tax) {
            ch_mag_input = EPITAXMAG_TAX.out.filtered
        } else {
            // Reuse filtered reads from a previous run
            def rn = params.run_name ?: file(params.input).name
            def od = params.outdir ?: "${projectDir}/results/${rn}"
            def filt_dir = params.filtered_dir ?: "${od}/04_filter_chopper"

            log.info "  MAG: reading filtered reads from ${filt_dir}"

            ch_mag_input = Channel
                .fromPath("${filt_dir}/*.filtered.fastq.gz", checkIfExists: true)
                .map { f -> tuple(f.baseName.replaceAll(/\.filtered\.fastq$/, ''), f) }
        }

        EPITAXMAG_MAG(ch_mag_input)
    }
}


// Emit a per-run params.json next to the reports so the HTML can render
// a fully reproducible "Run parameters" panel (which classifier, which
// DB versions, which thresholds, etc.). Runs once at the end regardless
// of whether the pipeline completed cleanly. DSL2 syntax: no `=`.
workflow.onComplete {
    try {
        def rn = params.run_name ?: file(params.input).name
        def od = params.outdir ?: "${projectDir}/results/${rn}"
        def reports_dir = file("${od}/00_reports")
        reports_dir.mkdirs()

        def safe = { v ->
            if (v == null) return null
            if (v instanceof Boolean || v instanceof Number) return v
            return v.toString()
        }

        def payload = [
            run: [
                name        : rn,
                started_at  : workflow.start?.toString(),
                completed_at: workflow.complete?.toString(),
                duration    : workflow.duration?.toString(),
                success     : workflow.success,
                exit_status : workflow.exitStatus,
                command_line: workflow.commandLine,
                nextflow    : workflow.nextflow?.version?.toString(),
                project_dir : workflow.projectDir?.toString(),
                launch_dir  : workflow.launchDir?.toString(),
                work_dir    : workflow.workDir?.toString(),
                profile     : workflow.profile,
                revision    : workflow.revision,
                session_id  : workflow.sessionId?.toString(),
            ],
            params: params.collectEntries { k, v -> [(k.toString()): safe(v)] }
        ]

        def json = groovy.json.JsonOutput.prettyPrint(
            groovy.json.JsonOutput.toJson(payload))
        file("${reports_dir}/params.json").text = json
        log.info "Wrote run parameters to ${reports_dir}/params.json"
    } catch (Throwable t) {
        log.warn "Could not write params.json: ${t.message}"
    }
}
