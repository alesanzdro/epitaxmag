#!/usr/bin/env nextflow

/*
 * ══════════════════════════════════════════════════════════════════
 *  EpiTaxMAG — Pipeline Metagenómica Nanopore
 *  Fase TAX (EpiTax): QC · Trimming · Filtrado · Perfilado taxonómico
 *  Fase MAG (EpiTaxMAG): Ensamblaje · Binning · Anotación · AMR
 *  EPIMOL — FISABIO
 * ══════════════════════════════════════════════════════════════════
 *
 *  USO:
 *    # Solo perfilado taxonómico (EpiTax)
 *    nextflow run main.nf --input data/260225_EPIM230 -profile gru
 *
 *    # Flujo completo: taxonómico + ensamblaje MAGs (EpiTaxMAG)
 *    nextflow run main.nf --input data/260225_EPIM230 -profile gru --run_assembly
 *
 *    # Solo ensamblaje MAGs (requiere reads filtrados previos)
 *    nextflow run main.nf --input data/260225_EPIM230 -profile gru --run_tax false --run_assembly
 *
 *    # Resume tras interrupción
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
// AYUDA
// ══════════════════════════════════════════════════════════════════

if (params.help) {
    log.info """
    EpiTaxMAG — Pipeline Metagenómica Nanopore
    EPIMOL — FISABIO
    ──────────────────────────────────────────

    USO:
      nextflow run main.nf --input <dir_fastqs> -profile <profile> [opciones]

    OPCIONES PRINCIPALES:
      --input            Directorio con ficheros .fastq.gz (REQUERIDO)
      --db_root          Raíz de bases de datos [default: ${projectDir}/databases]
      --outdir           Directorio de salida [default: results/<run_name>]
      --run_name         Nombre del run [default: nombre del dir de input]

    FASES:
      --run_tax          Ejecutar fase TAX [default: true]
      --run_assembly     Ejecutar fase MAG [default: false]

    FILTRADO:
      --min_quality      Calidad mínima Chopper [default: 10]
      --min_length       Longitud mínima (bp) [default: 1000]
      --host_genome      Genoma del huésped para HOST_REMOVAL (.fasta)
                         Si no se indica, el paso se omite

    TAXONOMÍA:
      --kraken2_confidence  Confianza Kraken2 [default: 0.1]
      --genome_size_mb      Tamaño genoma promedio Mb [default: 4.5]
      --min_coverage_mag    Cobertura mínima MAGs [default: 30]

    VERIFICACIÓN FENOTÍPICA (minimap2 vs type-strains):
      --phenotypic_targets  TSV con organismos objetivo + accesiones NCBI
                            (columnas: organism, accession) [default: activado]
                            Editar assets/phenotypic_targets/target_organisms.tsv
                            para añadir/quitar organismos.
                            Para desactivar: --phenotypic_targets false

    FASTQ SCREEN:
      --fastqscreen_conf    Config de FastQ Screen (.conf)

    BASES DE DATOS:
      --db_root             Raiz de bases de datos [default: databases/]
      --skip_db_setup       No descargar DBs automaticamente [default: false]

    OPCIONES MAG:
      --skip_gtdbtk         Saltar GTDB-Tk (equipos con <64 GB RAM) [default: false]
      --assemblies_dir      Reutilizar ensamblajes existentes (salta MetaFlye)
      --filtered_dir        Reutilizar reads filtrados (salta TAX)

    PROFILES (obligatorio elegir uno):
      -profile gru          Servidor GRU: 48 threads, 200 GB RAM
      -profile pc           Desktop: 20 threads, 56 GB RAM
      -profile censalud     CenSalud: 24 threads, 60 GB RAM (i7-14700)
      -profile singularity  Usar contenedores Singularity (combinable)

    EJEMPLOS:
      # Perfilado taxonómico con conda (por defecto)
      nextflow run main.nf --input data/EPIM230 -profile censalud

      # Con Singularity
      nextflow run main.nf --input data/EPIM230 -profile censalud,singularity

      # Con eliminación de host humano
      nextflow run main.nf --input data/EPIM230 -profile gru \\
          --host_genome /path/to/GRCh38.fna

      # Especificar bases de datos
      nextflow run main.nf --input data/EPIM230 -profile censalud \\
          --db_root /disco/databases

      # Sin verificación fenotípica (desactivar minimap2 type-strains)
      nextflow run main.nf --input data/EPIM230 -profile gru \\
          --phenotypic_targets false

      # Resume tras interrupción
      nextflow run main.nf --input data/EPIM230 -profile censalud -resume

    SETUP INICIAL (bases de datos):
      bash bin/setup_databases.sh --db-root /ruta/databases
      bash bin/setup_databases.sh --db-root /ruta/databases --check
    """.stripIndent()
    System.exit(0)
}


// ══════════════════════════════════════════════════════════════════
// VALIDACIÓN
// ══════════════════════════════════════════════════════════════════

if (!params.input) {
    error "ERROR: Falta --input. Uso: nextflow run main.nf --input /ruta/a/fastqs -profile gru|pc"
}

def run_name    = params.run_name ?: file(params.input).name
def pipeline_id = params.run_assembly ? 'EpiTaxMAG' : 'EpiTax'

// ── Validar bases de datos TAX ──────────────────────────────────
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
if (db_missing) {
    log.warn """
    ┌─────────────────────────────────────────────────────────────┐
    │  BASES DE DATOS NO ENCONTRADAS                              │
    └─────────────────────────────────────────────────────────────┘
    ${db_missing.join('\n    ')}

    Solución: ejecutar el setup de bases de datos:
      bash bin/setup_databases.sh --db-root ${params.db_root}

    O especificar una ruta existente:
      nextflow run main.nf --db_root /ruta/a/databases ...
    """.stripIndent()
}

log.info """
  =========================================
   E p i T a x M A G   v1.0
  =========================================
   Nanopore Metagenomics Pipeline
   FISABIO - EPIMOL
  -----------------------------------------

  Run:       ${run_name}
  Input:     ${params.input}
  Output:    ${params.outdir ?: "${projectDir}/results/${run_name}"}
  DB Root:   ${params.db_root}
  Fase TAX:  ${params.run_tax ? 'SI' : 'NO'}
  Fase MAG:  ${params.run_assembly ? 'SI' : 'NO'}
  Filtros:   Q>=${params.min_quality}, len>=${params.min_length}bp
  MAG cov:   >=${params.min_coverage_mag}x (genoma ${params.genome_size_mb} Mb)
""".stripIndent()


// ══════════════════════════════════════════════════════════════════
// CANAL DE ENTRADA — descubrir muestras (excluye unclassified)
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
// WORKFLOW PRINCIPAL
// ══════════════════════════════════════════════════════════════════

workflow {

    // ── Setup: Descarga automatica de bases de datos ─────────────
    if (!params.skip_db_setup) {
        EPITAXMAG_SETUP()
    }

    // ── Fase TAX: Perfilado taxonómico (EpiTax) ──────────────────
    if (params.run_tax) {
        EPITAXMAG_TAX(ch_raw_fastq)
    }

    // ── Fase MAG: Ensamblaje + MAGs (EpiTaxMAG) ─────────────────
    if (params.run_assembly) {
        if (params.run_tax) {
            ch_mag_input = EPITAXMAG_TAX.out.filtered
        } else {
            // Buscar reads filtrados de un run previo
            def rn = params.run_name ?: file(params.input).name
            def od = params.outdir ?: "${projectDir}/results/${rn}"
            def filt_dir = params.filtered_dir ?: "${od}/04_filter_chopper"

            log.info "  MAG: buscando reads filtrados en ${filt_dir}"

            ch_mag_input = Channel
                .fromPath("${filt_dir}/*.filtered.fastq.gz", checkIfExists: true)
                .map { f -> tuple(f.baseName.replaceAll(/\.filtered\.fastq$/, ''), f) }
        }

        EPITAXMAG_MAG(ch_mag_input)
    }
}
