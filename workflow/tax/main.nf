#!/usr/bin/env nextflow

/*
 * ══════════════════════════════════════════════════════════════════
 *  EpiTax — Subworkflow de Perfilado Taxonómico
 *  QC · Trimming · Filtrado · Host Removal · FastQ Screen ·
 *  Taxonomía · AMR screening · Reporte
 *  EPIMOL — FISABIO
 * ══════════════════════════════════════════════════════════════════
 *
 *  PIPELINE (14 pasos):
 *   01. QC raw (NanoPlot + FastQC + CHECK_SQCORE)
 *   02. FastQ Screen (sobre reads RAW — diagnóstico de vectores, host, etc.)
 *   03. Porechop ABI (trimming)
 *   04. Chopper (filtrado Q + longitud)
 *   05. HOST_REMOVAL (opcional, --host_genome)
 *   06. QC filtrado (NanoPlot + FastQC + CHECK_SQCORE)
 *   07. Kraken2
 *   08. Bracken
 *   09. Kaiju
 *   10. Sylph
 *   11. KMA AMR screening (ResFinder desde reads)
 *   12. Verificación fenotípica (minimap2 vs type-strains)
 *   13. MultiQC
 *   14. Reporte HTML EpiTax (Plotly, incluye screening + fenotípico)
 */


// ── Variables derivadas de params ──────────────────────────────────
def run_name = params.run_name ?: (params.input ? file(params.input).name : 'unknown')
def outdir   = params.outdir   ?: "${projectDir}/results/${run_name}"


// ══════════════════════════════════════════════════════════════════
// PASO 01: QC reads crudos
// ══════════════════════════════════════════════════════════════════

process NANOPLOT_RAW {
    tag "${sample}"
    publishDir "${outdir}/01_qc_raw/nanoplot/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    path("*"), emit: nanoplot_files

    script:
    """
    NanoPlot --fastq ${fastq} \\
        --outdir . \\
        --prefix ${sample}_ \\
        --threads ${task.cpus} \\
        --loglength --N50 \\
        --title "${sample} - Raw"
    """
}

process FASTQC_RAW {
    tag "${sample}"
    publishDir "${outdir}/01_qc_raw/fastqc", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    path("*"), emit: fastqc_files

    script:
    """
    fastqc ${fastq} --outdir . --threads ${task.cpus} --quiet
    """
}

process CHECK_SQCORE_RAW {
    publishDir "${outdir}/01_qc_raw", mode: 'copy'

    input:
    path(fastq_dir)

    output:
    path("qc_raw_metrics.csv"),   emit: csv
    path("*_report.pdf"),         emit: pdf, optional: true

    script:
    """
    python3 ${params.scripts_dir}/01_check_sqcore_v2.py ${fastq_dir} -t ${task.cpus}
    mv ${fastq_dir}/analysis_with_q_custom.csv qc_raw_metrics.csv || true
    python3 ${params.scripts_dir}/generate_qc_report.py qc_raw_metrics.csv ${run_name}_qc_raw_report.pdf || true
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 02: FASTQ SCREEN — screening contaminación (sobre reads RAW)
// Corre sobre rawdata para detectar vectores, adaptadores, host,
// controles y patógenos ANTES de que Porechop los elimine.
// ══════════════════════════════════════════════════════════════════

process FASTQSCREEN {
    tag "${sample}"
    publishDir "${outdir}/02_fastqscreen_raw", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    path("*_screen.txt"), emit: screen_txt
    path("*"),            emit: screen_files

    script:
    """
    # FastQ Screen: el _screen.txt se genera antes del HTML.
    # Si falla la generacion del HTML (template missing en container),
    # toleramos el error si el .txt existe.
    fastq_screen --conf ${params.fastqscreen_conf} \\
        --outdir . \\
        --threads ${task.cpus} \\
        --aligner minimap2 \\
        --minimap2 '-x map-ont' \\
        ${fastq} || true

    # Verificar que el output esencial existe
    ls *_screen.txt > /dev/null 2>&1 || { echo "ERROR: no _screen.txt generated"; exit 1; }
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 03: PORECHOP_ABI — trimming adaptadores
// ══════════════════════════════════════════════════════════════════

process PORECHOP {
    tag "${sample}"
    publishDir "${outdir}/03_trim_porechop", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.trimmed.fastq.gz"), emit: trimmed

    script:
    """
    porechop_abi \\
        -i ${fastq} \\
        -o ${sample}.trimmed.fastq.gz \\
        -t ${task.cpus}
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 04: CHOPPER — filtrado calidad + longitud
// ══════════════════════════════════════════════════════════════════

process CHOPPER {
    tag "${sample}"
    publishDir "${outdir}/04_filter_chopper", mode: 'copy'

    input:
    tuple val(sample), path(trimmed)

    output:
    tuple val(sample), path("${sample}.filtered.fastq.gz"), emit: filtered

    script:
    """
    gunzip -c ${trimmed} \\
        | chopper --quality ${params.min_quality} \\
                  --minlength ${params.min_length} \\
                  --threads ${task.cpus} \\
        | gzip > ${sample}.filtered.fastq.gz
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 05: HOST_REMOVAL — depleción de reads del hospedador
// Solo se ejecuta si params.host_genome != null
// ══════════════════════════════════════════════════════════════════

process HOST_REMOVAL {
    tag "${sample}"
    publishDir "${outdir}/05_host_removal", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.clean.fastq.gz"), emit: clean
    path("${sample}.host_stats.txt"), emit: stats

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${params.host_genome} ${fastq} \\
        2>/dev/null \\
        | samtools view -bSu -f 4 -@ 4 - \\
        | samtools sort -n -@ 4 - \\
        | samtools fastq -@ 2 - \\
        | gzip > ${sample}.clean.fastq.gz

    total=\$(zcat ${fastq} | awk 'NR%4==1' | wc -l)
    clean=\$(zcat ${sample}.clean.fastq.gz | awk 'NR%4==1' | wc -l)
    removed=\$((total - clean))
    pct=\$(echo "scale=1; \${clean}*100/\${total}" | bc)
    echo -e "${sample}\\t\${total}\\t\${removed}\\t\${clean}\\t\${pct}%" \\
        > ${sample}.host_stats.txt
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 06: QC reads filtrados
// ══════════════════════════════════════════════════════════════════

process NANOPLOT_FILTERED {
    tag "${sample}"
    publishDir "${outdir}/06_qc_filtered/nanoplot/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    path("*"), emit: nanoplot_files

    script:
    """
    NanoPlot --fastq ${fastq} \\
        --outdir . \\
        --prefix ${sample}_ \\
        --threads ${task.cpus} \\
        --loglength --N50 \\
        --title "${sample} - Filtered"
    """
}

process FASTQC_FILTERED {
    tag "${sample}"
    publishDir "${outdir}/06_qc_filtered/fastqc", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    path("*"), emit: fastqc_files

    script:
    """
    fastqc ${fastq} --outdir . --threads ${task.cpus} --quiet
    """
}

process CHECK_SQCORE_FILTERED {
    publishDir "${outdir}/06_qc_filtered", mode: 'copy'

    input:
    path(filtered_fastqs)

    output:
    path("qc_filtered_metrics.csv"),  emit: csv
    path("*_report.pdf"),             emit: pdf, optional: true

    script:
    def tag = params.host_genome ? ".clean" : ".filtered"
    """
    mkdir -p fastqs
    cp ${filtered_fastqs} fastqs/
    python3 ${params.scripts_dir}/01_check_sqcore_v2.py fastqs/ -t ${task.cpus} --tagdelete "${tag}"
    mv fastqs/analysis_with_q_custom.csv qc_filtered_metrics.csv || true
    python3 ${params.scripts_dir}/generate_qc_report.py qc_filtered_metrics.csv ${run_name}_qc_filtered_report.pdf || true
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 07: KRAKEN2 — clasificación k-mers
// ══════════════════════════════════════════════════════════════════

process KRAKEN2 {
    tag "${sample}"
    cpus params.threads_heavy ?: 20
    publishDir "${outdir}/07_tax_kraken2", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.kraken2.report"), emit: report
    path("${sample}.kraken2.out"),                       emit: output

    script:
    def mmap = (!task.memory || task.memory.toGiga() < 100) ? '--memory-mapping' : ''
    """
    kraken2 \\
        --db ${params.db_root}/k2_pluspf_20251015 \\
        ${mmap} \\
        --threads ${task.cpus} \\
        --confidence ${params.kraken2_confidence} \\
        --output ${sample}.kraken2.out \\
        --report ${sample}.kraken2.report \\
        ${fastq}
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 08: BRACKEN — estimación abundancias
// ══════════════════════════════════════════════════════════════════

process BRACKEN {
    tag "${sample}"
    cpus 2
    memory '4 GB'
    publishDir "${outdir}/08_tax_bracken", mode: 'copy'

    input:
    tuple val(sample), path(kraken_report)

    output:
    path("${sample}.bracken.*"), emit: bracken_files

    script:
    """
    for level in S G; do
        bracken \\
            -d ${params.db_root}/k2_pluspf_20251015 \\
            -i ${kraken_report} \\
            -o ${sample}.bracken.\${level}.txt \\
            -w ${sample}.bracken.\${level}.report \\
            -r 150 -l \${level} -t 10 || true
    done
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 09: KAIJU — clasificación proteica
// ══════════════════════════════════════════════════════════════════

process KAIJU {
    tag "${sample}"
    cpus params.threads_heavy ?: 20
    publishDir "${outdir}/09_tax_kaiju", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    path("${sample}.kaiju.*"), emit: kaiju_files

    script:
    """
    kaiju \\
        -t ${params.db_root}/kaiju/nodes.dmp \\
        -f ${params.db_root}/kaiju/kaiju_db_refseq_ref.fmi \\
        -i ${fastq} \\
        -o ${sample}.kaiju.out \\
        -z ${task.cpus}

    kaiju2table -t ${params.db_root}/kaiju/nodes.dmp -n ${params.db_root}/kaiju/names.dmp \\
        -r species -o ${sample}.kaiju.species.tsv ${sample}.kaiju.out

    kaiju2table -t ${params.db_root}/kaiju/nodes.dmp -n ${params.db_root}/kaiju/names.dmp \\
        -r genus -o ${sample}.kaiju.genus.tsv ${sample}.kaiju.out
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 10: SYLPH — perfilado ANI (GTDB + Fungi + Viral)
// ══════════════════════════════════════════════════════════════════

process SYLPH {
    cpus params.threads_heavy ?: 20
    publishDir "${outdir}/10_tax_sylph", mode: 'copy'

    input:
    path(filtered_fastqs)

    output:
    path("sylph_profile_all.tsv"),   emit: profile
    path("*.sylph.tsv"),             emit: per_sample, optional: true
    path("*_report.pdf"),            emit: pdf, optional: true

    script:
    def ext = params.host_genome ? ".clean.fastq.gz" : ".filtered.fastq.gz"
    """
    sylph profile \\
        ${params.db_root}/sylph/gtdb-r226-c200-dbv1.syldb ${params.db_root}/sylph/fungi-refseq-2024-07-25-c200-v0.3.syldb ${params.db_root}/sylph/imgvr_c200_v0.3.0.syldb \\
        ${filtered_fastqs} \\
        -t ${task.cpus} \\
        > sylph_profile_all.tsv

    head -1 sylph_profile_all.tsv > _header.tsv
    for fq in ${filtered_fastqs}; do
        sample=\$(basename \$fq ${ext})
        cp _header.tsv \${sample}.sylph.tsv
        grep "\${sample}" sylph_profile_all.tsv >> \${sample}.sylph.tsv || true
    done
    rm _header.tsv

    python3 ${params.scripts_dir}/sylph_metagenomics_report.py \\
        sylph_profile_all.tsv ${run_name}_sylph_report.pdf || true
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 11: KMA_READS — AMR screening desde reads (ResFinder)
// Detección rápida de genes de resistencia, virulencia, plásmidos
// sin necesidad de ensamblaje
// ══════════════════════════════════════════════════════════════════

process KMA_READS {
    tag "${sample}"
    conda "${projectDir}/envs/nf-kma.yml"
    publishDir "${outdir}/11_amr_kma", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    path("${sample}.res"),  emit: results
    path("${sample}.*"),    emit: all_files

    script:
    """
    # Indexar la DB ResFinder combinada en el workdir (3 MB, < 1 seg)
    kma index -i ${params.db_root}/kma_resfinder/resfinder_db/all.fsa -o resfinder_idx 2>/dev/null

    # Mapeo de reads contra genes de resistencia
    kma \\
        -i ${fastq} \\
        -o ${sample} \\
        -t_db resfinder_idx \\
        -ont \\
        -t ${task.cpus} \\
        -1t1 \\
        -mem_mode \\
        2>/dev/null || true

    # Asegurar que el fichero .res existe aunque no haya hits
    touch ${sample}.res
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 12: VERIFICACIÓN FENOTÍPICA — minimap2 vs type-strains
// Opcional: solo si params.phenotypic_targets está definido.
// Descarga genomas de referencia (type-strain) de NCBI y mapea
// los reads filtrados para verificar presencia/ausencia.
// Reporta breadth of coverage a ≥1x, ≥10x, ≥30x de profundidad.
// ══════════════════════════════════════════════════════════════════

process BUILD_PHENOTYPIC_DB {
    publishDir "${outdir}/12_phenotypic_minimap2/db", mode: 'copy'

    input:
    path(targets_tsv)

    output:
    path("phenotypic_targets.fa"),  emit: fasta
    path("phenotypic_targets.mmi"), emit: index

    script:
    """
    python3 ${params.scripts_dir}/build_phenotypic_db.py \\
        --targets ${targets_tsv} \\
        --output phenotypic_targets.fa

    minimap2 -d phenotypic_targets.mmi phenotypic_targets.fa
    """
}


process MINIMAP2_PHENOTYPIC {
    tag "${sample}"
    publishDir "${outdir}/12_phenotypic_minimap2", mode: 'copy', pattern: "*_phenotypic.tsv"

    input:
    tuple val(sample), path(fastq)
    path(fasta)
    path(index)

    output:
    path("${sample}_phenotypic.tsv"), emit: coverage

    script:
    """
    # Map reads against type-strain genomes
    minimap2 -ax map-ont -t ${task.cpus} ${index} ${fastq} \\
        | samtools sort -@ 2 -o ${sample}.bam
    samtools index ${sample}.bam

    # Per-contig depth-stratified coverage (streaming awk — single pass)
    samtools depth -a ${sample}.bam | awk -F'\\t' '{
        ref=\$1; depth=\$3
        total[ref]++
        if (depth >= 1) d1[ref]++
        if (depth >= 10) d10[ref]++
        if (depth >= 30) d30[ref]++
        dsum[ref] += depth
    } END {
        print "reference\\tgenome_size\\tbases_1x\\tbases_10x\\tbases_30x\\tdepth_sum"
        for (r in total) {
            printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n", r, total[r], d1[r]+0, d10[r]+0, d30[r]+0, dsum[r]
        }
    }' > depth_stats.tsv

    # Mapped reads per contig
    samtools idxstats ${sample}.bam > idxstats.tsv

    # Aggregate by species (chromosome + plasmids → one row per species)
    python3 ${params.scripts_dir}/aggregate_phenotypic_coverage.py \\
        --depth-stats depth_stats.tsv \\
        --idxstats idxstats.tsv \\
        --sample ${sample} \\
        --output ${sample}_phenotypic.tsv
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 13: MULTIQC — informe integrado
// ══════════════════════════════════════════════════════════════════

process MULTIQC {
    cpus 2
    memory '4 GB'
    publishDir "${outdir}/13_multiqc", mode: 'copy'

    input:
    val(ready)

    output:
    path("${run_name}_multiqc*"), emit: report

    script:
    """
    multiqc ${outdir} \\
        --outdir . \\
        --title "EpiTax — ${run_name}" \\
        --filename ${run_name}_multiqc_report \\
        --config ${projectDir}/multiqc_config.yaml \\
        --force \\
        --quiet
    """
}


// ══════════════════════════════════════════════════════════════════
// PASO 14: REPORTE FINAL HTML EpiTax (último paso)
// ══════════════════════════════════════════════════════════════════

process FINAL_REPORT {
    cpus 2
    memory '4 GB'
    publishDir "${outdir}/14_epitax_report", mode: 'copy'

    input:
    path(qc_raw_csv)
    path(qc_filtered_csv)
    path(taxonomy_files)
    path(sylph_profile)
    path(fastqscreen_files)
    path(kma_files)
    path(phenotypic_files)

    output:
    path("${run_name}_final_report.html"), emit: report

    script:
    """
    # Combine per-sample phenotypic TSVs if present
    PHENO_ARG=""
    if ls *_phenotypic.tsv 1>/dev/null 2>&1; then
        head -1 \$(ls *_phenotypic.tsv | head -1) > phenotypic_combined.tsv
        for f in *_phenotypic.tsv; do tail -n +2 "\$f" >> phenotypic_combined.tsv; done
        PHENO_ARG="--phenotypic phenotypic_combined.tsv"
    fi

    python3 ${params.scripts_dir}/generate_html_summary.py \\
        --qc-raw ${qc_raw_csv} \\
        --qc-filtered ${qc_filtered_csv} \\
        --sylph ${sylph_profile} \\
        --run-name ${run_name} \\
        --output ${run_name}_final_report.html \\
        --genome-size ${params.genome_size_mb} \\
        --min-coverage ${params.min_coverage_mag} \\
        --min-quality ${params.min_quality} \\
        --min-length ${params.min_length} \\
        \$PHENO_ARG
    """
}


// ══════════════════════════════════════════════════════════════════
// WORKFLOW TAX — Perfilado Taxonómico Completo (EpiTax)
// ══════════════════════════════════════════════════════════════════

workflow TAX {

    take:
    ch_raw_fastq    // tuple(sample, fastq)

    main:

    // ── Paso 01: QC raw ────────────────────────────────────────
    NANOPLOT_RAW(ch_raw_fastq)
    FASTQC_RAW(ch_raw_fastq)
    CHECK_SQCORE_RAW(Channel.of(file(params.input)))

    // ── Paso 02: FastQ Screen (sobre reads RAW) ───────────────
    // Diagnóstico de vectores, host, controles, patógenos
    // ANTES de que Porechop elimine adaptadores/vectores
    FASTQSCREEN(ch_raw_fastq)

    // ── Paso 03: Porechop — trimming adaptadores ───────────────
    PORECHOP(ch_raw_fastq)

    // ── Paso 04: Chopper — filtrado Q + longitud ───────────────
    CHOPPER(PORECHOP.out.trimmed)

    // ── Paso 05: Host removal (opcional) ──────────────────────
    if (params.host_genome) {
        HOST_REMOVAL(CHOPPER.out.filtered)
        ch_clean = HOST_REMOVAL.out.clean
    } else {
        ch_clean = CHOPPER.out.filtered
    }

    // ── Paso 06: QC filtrados (sobre reads de Chopper) ────────
    NANOPLOT_FILTERED(CHOPPER.out.filtered)
    FASTQC_FILTERED(CHOPPER.out.filtered)

    ch_all_clean = ch_clean
        .map { sample, fastq -> fastq }
        .collect()

    CHECK_SQCORE_FILTERED(ch_all_clean)

    // ── Paso 07-08: Kraken2 → Bracken (sobre reads limpios) ──
    KRAKEN2(ch_clean)
    BRACKEN(KRAKEN2.out.report)

    // ── Paso 09: Kaiju (sobre reads limpios) ─────────────────
    KAIJU(ch_clean)

    // ── Paso 10: Sylph (sobre reads limpios) ─────────────────
    SYLPH(ch_all_clean)

    // ── Paso 11: KMA AMR screening (sobre reads limpios) ─────
    KMA_READS(ch_clean)

    // ── Paso 12: Verificación fenotípica (minimap2 vs type-strains)
    // Descarga genomas de referencia (NCBI) y mapea reads para
    // verificar presencia/ausencia de organismos objetivo.
    // Activado por defecto. Desactivar: --phenotypic_targets false
    ch_phenotypic_trigger = Channel.empty()
    if (params.phenotypic_targets && params.phenotypic_targets != 'false') {
        BUILD_PHENOTYPIC_DB(Channel.of(file(params.phenotypic_targets)))
        MINIMAP2_PHENOTYPIC(
            ch_clean,
            BUILD_PHENOTYPIC_DB.out.fasta.first(),
            BUILD_PHENOTYPIC_DB.out.index.first()
        )
        ch_phenotypic = MINIMAP2_PHENOTYPIC.out.coverage.collect()
        ch_phenotypic_trigger = ch_phenotypic
    } else {
        ch_phenotypic = Channel.value([])
    }

    // ── Paso 13: MultiQC — esperar a que TODO termine ─────────
    Channel.empty()
        .mix(
            FASTQC_RAW.out.fastqc_files.collect(),
            FASTQC_FILTERED.out.fastqc_files.collect(),
            FASTQSCREEN.out.screen_files.collect(),
            BRACKEN.out.bracken_files.collect(),
            KAIJU.out.kaiju_files.collect(),
            SYLPH.out.profile,
            KMA_READS.out.all_files.collect(),
            ch_phenotypic_trigger
        )
        .collect()
        .map { 'ready' }
        .set { ch_multiqc_trigger }

    MULTIQC(ch_multiqc_trigger)

    // ── Paso 14: Reporte final HTML EpiTax ─────────────────
    ch_taxonomy = Channel.empty()
        .mix(
            KRAKEN2.out.report.map { sample, report -> report },
            BRACKEN.out.bracken_files.flatten(),
            KAIJU.out.kaiju_files.flatten()
        )
        .collect()

    ch_fastqscreen = FASTQSCREEN.out.screen_txt.collect()
    ch_kma = KMA_READS.out.results.collect()

    FINAL_REPORT(
        CHECK_SQCORE_RAW.out.csv,
        CHECK_SQCORE_FILTERED.out.csv,
        ch_taxonomy,
        SYLPH.out.profile,
        ch_fastqscreen,
        ch_kma,
        ch_phenotypic
    )

    emit:
    filtered        = ch_clean
    qc_raw_csv      = CHECK_SQCORE_RAW.out.csv
    qc_filtered_csv = CHECK_SQCORE_FILTERED.out.csv
    sylph_profile   = SYLPH.out.profile
    report          = FINAL_REPORT.out.report
}
