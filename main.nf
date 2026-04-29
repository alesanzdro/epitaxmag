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
      --skip_gtdbtk         Skip GTDB-Tk (machines with <50 GB RAM, uses Sourmash)
      --assemblies_dir      Reuse existing assemblies (skips MetaFlye)
      --filtered_dir        Reuse pre-filtered reads (skips TAX)

    PROFILES (one is required):
      -profile gru          Server: 48 threads, 200 GB RAM
      -profile pc           Desktop: 20 threads, 56 GB RAM
      -profile censalud     CenSalud: 24 threads, 60 GB RAM (i7-14700)
      -profile singularity  Use Singularity containers (combinable)

    EXAMPLES:
      # Taxonomic profiling with conda (default)
      nextflow run main.nf --input data/EPIM230 -profile censalud

      # With Singularity
      nextflow run main.nf --input data/EPIM230 -profile censalud,singularity

      # With human host removal
      nextflow run main.nf --input data/EPIM230 -profile gru \\
          --host_genome /path/to/GRCh38.fna

      # Custom database root
      nextflow run main.nf --input data/EPIM230 -profile censalud \\
          --db_root /disk/databases

      # Disable phenotypic verification
      nextflow run main.nf --input data/EPIM230 -profile gru \\
          --phenotypic_targets false

      # Resume after interruption
      nextflow run main.nf --input data/EPIM230 -profile censalud -resume
    """.stripIndent()
    System.exit(0)
}


// ══════════════════════════════════════════════════════════════════
// VALIDATION
// ══════════════════════════════════════════════════════════════════

if (!params.input) {
    error "ERROR: --input is required. Usage: nextflow run main.nf --input /path/to/fastqs -profile <gru|pc|censalud>"
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
