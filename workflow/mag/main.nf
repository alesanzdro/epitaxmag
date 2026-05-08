#!/usr/bin/env nextflow

/*
 * EpiTaxMAG — MAG core subworkflow
 * MetaFlye + Medaka + QUAST + Binning (MetaBAT2 + MaxBin2 + SemiBin2)
 * + DAS Tool + CheckM2 + GTDB-Tk/Sourmash + AMRFinderPlus + geNomad
 * + ABRicate + MOB-suite + IntegronFinder + Bakta
 * EPIMOL — FISABIO
 *
 * PIPELINE (17 steps):
 *  01. MetaFlye  — de novo metagenomic assembly
 *  02. Medaka    — polishing (SUP v5 model)
 *  03. QUAST     — assembly statistics
 *  04. MAP_READS — read-to-assembly mapping (coverage)
 *  05. MetaBAT2  — composition + coverage binning
 *  06. MaxBin2   — EM-based binning
 *  07. SemiBin2  — deep-learning binning (long reads)
 *  08. DAS Tool  — bin refinement / consensus
 *  09. CheckM2   — MAG completeness and contamination
 *  10. GTDB-Tk   — taxonomic classification of MAGs (or SOURMASH_CLASSIFY)
 *  11. AMRFinderPlus — AMR + virulence + stress in MAGs
 *  12. geNomad   — plasmid + virus detection on assemblies
 *  13. ABRicate  — virulence (VFDB) + AMR (CARD) + PlasmidFinder
 *  14. MOB-suite — plasmid typing
 *  15. IntegronFinder — integron detection per MAG
 *  16. Bakta     — functional annotation per MAG
 *  17. AMR-pathogen integration + final HTML/Excel report
 */

def run_name = params.run_name ?: (params.input ? file(params.input).name : 'unknown')
def outdir   = params.outdir   ?: "${projectDir}/results/${run_name}"


// ======================================================================
// STEP 01: METAFLYE — de novo metagenomic assembly
// ======================================================================

process METAFLYE {
    tag "${sample}"
    publishDir "${outdir}/15_assembly_flye/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.assembly.fasta"), emit: assembly
    path("${sample}.assembly_info.txt"),                 emit: info

    script:
    """
    flye --meta ${params.flye_mode} ${fastq} \
        --out-dir flye_out \
        --threads ${task.cpus}

    cp flye_out/assembly.fasta ${sample}.assembly.fasta
    cp flye_out/assembly_info.txt ${sample}.assembly_info.txt
    """
}


// ======================================================================
// STEP 02: MEDAKA — assembly polishing
// ======================================================================

process MEDAKA {
    tag "${sample}"
    publishDir "${outdir}/16_polish_medaka", mode: 'copy'

    input:
    tuple val(sample), path(assembly), path(fastq)

    output:
    tuple val(sample), path("${sample}.polished.fasta"), emit: polished

    script:
    """
    # 1. Align reads to the assembly
    minimap2 -ax map-ont -t ${task.cpus} ${assembly} ${fastq} \
        | samtools sort -@ 4 -o aligned.bam
    samtools index aligned.bam

    # 2. Run Medaka neural-network inference
    medaka inference aligned.bam medaka_out \
        --model ${params.medaka_model} \
        --bam_workers 2

    # 3. Build the polished consensus
    medaka sequence medaka_out ${assembly} ${sample}.polished.fasta

    # Cleanup
    rm -f aligned.bam aligned.bam.bai
    """
}


// ======================================================================
// STEP 03: QUAST — assembly statistics
// ======================================================================

process QUAST {
    publishDir "${outdir}/17_qc_assembly", mode: 'copy'

    input:
    path(assemblies)
    val(labels)

    output:
    path("quast_results/"), emit: report

    script:
    def label_str = labels.join(',')
    """
    quast ${assemblies} \
        --labels "${label_str}" \
        --output-dir quast_results \
        --threads ${task.cpus} \
        --min-contig 500
    """
}


// ======================================================================
// STEP 04: MAP_READS — map reads to assembly for coverage
// ======================================================================

process MAP_READS {
    tag "${sample}"
    publishDir "${outdir}/18_coverage", mode: 'copy'

    input:
    tuple val(sample), path(assembly), path(fastq)

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam
    tuple val(sample), path(assembly),                                                  emit: assembly

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${assembly} ${fastq} \
        | samtools sort -@ 4 -o ${sample}.sorted.bam
    samtools index ${sample}.sorted.bam
    """
}


// ======================================================================
// STEP 05: METABAT2 — composition + coverage binning
// ======================================================================

process METABAT2 {
    tag "${sample}"
    publishDir "${outdir}/19_binning_metabat2/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly), path(bam), path(bai)

    output:
    tuple val(sample), path("metabat_bins/"), emit: bins

    script:
    """
    jgi_summarize_bam_contig_depths \
        --outputDepth depth.txt \
        ${bam}

    mkdir -p metabat_bins
    metabat2 \
        -i ${assembly} \
        -a depth.txt \
        -o metabat_bins/${sample}_metabat \
        -t ${task.cpus} \
        --minContig 2500

    touch metabat_bins/.done
    """
}


// ======================================================================
// STEP 06: MAXBIN2 — expectation-maximization binning
// ======================================================================

process MAXBIN2 {
    tag "${sample}"
    publishDir "${outdir}/20_binning_maxbin2/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly), path(fastq)

    output:
    tuple val(sample), path("maxbin_bins/"), emit: bins

    script:
    """
    mkdir -p maxbin_bins
    run_MaxBin.pl \
        -contig ${assembly} \
        -reads ${fastq} \
        -out maxbin_bins/${sample}_maxbin \
        -thread ${task.cpus} \
        -min_contig_length 2500 || true

    # MaxBin2 genera .fasta — renombrar a .fa
    for f in maxbin_bins/*.fasta; do
        [ -f "\$f" ] && mv "\$f" "\${f%.fasta}.fa"
    done
    touch maxbin_bins/.done
    """
}


// ======================================================================
// STEP 07: SEMIBIN2 — deep-learning binning (long reads)
// ======================================================================

process SEMIBIN2 {
    tag "${sample}"
    publishDir "${outdir}/21_binning_semibin2/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly), path(bam), path(bai)

    output:
    tuple val(sample), path("semibin_bins/"), emit: bins

    script:
    """
    SemiBin2 single_easy_bin \
        -i ${assembly} \
        -b ${bam} \
        -o semibin_out \
        --sequencing-type long_read \
        --environment global \
        -t ${task.cpus}

    mkdir -p semibin_bins
    # SemiBin2 genera .fa.gz — descomprimir a .fa
    for gz in semibin_out/output_bins/*.fa.gz; do
        [ -f "\$gz" ] && gunzip -c "\$gz" > "semibin_bins/\$(basename "\$gz" .gz)"
    done
    # Por si genera .fa sin comprimir
    cp semibin_out/output_bins/*.fa semibin_bins/ 2>/dev/null || true
    touch semibin_bins/.done
    """
}


// ======================================================================
// STEP 08: DAS_TOOL — bin refinement / consensus
// ======================================================================

process DAS_TOOL {
    tag "${sample}"
    publishDir "${outdir}/22_binning_dastool/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly), path(metabat_bins), path(maxbin_bins), path(semibin_bins)

    output:
    tuple val(sample), path("${sample}_DASTool_bins/"), emit: refined_bins
    path("${sample}_DASTool_summary.tsv"),              emit: summary, optional: true

    script:
    """
    # Generate scaffolds2bin for each binner
    mkdir -p s2b

    # MetaBAT2 bins (Fasta_to_Contig2Bin genera 4 cols con metabat, cortar a 2)
    if ls ${metabat_bins}/*.fa 1>/dev/null 2>&1; then
        Fasta_to_Contig2Bin.sh -i ${metabat_bins} -e fa 2>/dev/null | cut -f1,4 > s2b/metabat2.tsv
    else
        touch s2b/metabat2.tsv
    fi

    # MaxBin2 bins
    if ls ${maxbin_bins}/*.fa 1>/dev/null 2>&1; then
        Fasta_to_Contig2Bin.sh -i ${maxbin_bins} -e fa > s2b/maxbin2.tsv 2>/dev/null
    else
        touch s2b/maxbin2.tsv
    fi

    # SemiBin2 bins
    if ls ${semibin_bins}/*.fa 1>/dev/null 2>&1; then
        Fasta_to_Contig2Bin.sh -i ${semibin_bins} -e fa > s2b/semibin2.tsv 2>/dev/null
    else
        touch s2b/semibin2.tsv
    fi

    # Preparar listas (solo binners con resultados)
    FILES=""
    LABELS=""
    for binner in metabat2 maxbin2 semibin2; do
        if [ -s "s2b/\${binner}.tsv" ]; then
            if [ -n "\$FILES" ]; then
                FILES="\${FILES},s2b/\${binner}.tsv"
                LABELS="\${LABELS},\${binner}"
            else
                FILES="s2b/\${binner}.tsv"
                LABELS="\${binner}"
            fi
        fi
    done

    mkdir -p ${sample}_DASTool_bins

    # Contar cuantos binners tienen resultados
    N_BINNERS=0
    [ -s s2b/metabat2.tsv ] && N_BINNERS=\$((N_BINNERS+1))
    [ -s s2b/semibin2.tsv ] && N_BINNERS=\$((N_BINNERS+1))

    if [ "\$N_BINNERS" -ge 2 ]; then
        # Dos o mas binners: DAS Tool refina
        DAS_Tool \
            -i "\$FILES" \
            -l "\$LABELS" \
            -c ${assembly} \
            -o ${sample} \
            --write_bins \
            --threads ${task.cpus} || true
    elif [ "\$N_BINNERS" -eq 1 ]; then
        # Solo un binner: copiar sus bins directamente (DAS Tool falla con 1 input)
        echo "Solo 1 binner con resultados, copiando bins directamente"
        for bdir in ${metabat_bins} ${maxbin_bins} ${semibin_bins}; do
            if ls \${bdir}/*.fa 1>/dev/null 2>&1; then
                cp \${bdir}/*.fa ${sample}_DASTool_bins/ 2>/dev/null || true
                break
            fi
        done
    else
        echo "Ningun binner produjo bins para ${sample}"
    fi

    touch ${sample}_DASTool_bins/.done
    """
}


// ======================================================================
// STEP 09: CHECKM2 — MAG completeness and contamination
// ======================================================================

process CHECKM2 {
    tag "${sample}"
    publishDir "${outdir}/23_binqc_checkm2/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bins_dir)

    output:
    // Emit a (sample, file) tuple so the workflow can join CHECKM2's
    // output with the other per-sample channels that feed
    // AMR_PATHOGEN_INTEGRATION. The filename `quality_report.tsv`
    // does not carry the sample id, so we keep it as the second
    // element of the tuple.
    tuple val(sample), path("quality_report.tsv"), emit: report, optional: true
    path("checkm2_out/"),                          emit: full_output, optional: true

    script:
    """
    # Contar bins .fa reales (no el .done marker)
    N_BINS=\$(ls ${bins_dir}/*.fa 2>/dev/null | wc -l)

    if [ "\$N_BINS" -gt 0 ]; then
        checkm2 predict \
            --input ${bins_dir} \
            --output-directory checkm2_out \
            --extension fa \
            --threads ${task.cpus} \
            --database_path ${params.db_root}/checkm2/CheckM2_database/uniref100.KO.1.dmnd

        cp checkm2_out/quality_report.tsv . 2>/dev/null || true
    else
        echo "No bins found for ${sample}, skipping CheckM2"
        mkdir -p checkm2_out
        echo -e "Name\tCompleteness\tContamination\tCompleteness_Model_Used" > quality_report.tsv
    fi
    """
}


// ======================================================================
// STEP 10: GTDB-Tk — taxonomic classification of MAGs
// ======================================================================

process GTDBTK {
    tag "${sample}"
    publishDir "${outdir}/24_taxonomy_gtdbtk/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bins_dir)

    output:
    path("gtdbtk_results/"),        emit: report
    path("${sample}_taxonomy.tsv"), emit: summary, optional: true

    script:
    """
    N_BINS=\$(ls ${bins_dir}/*.fa 2>/dev/null | wc -l)

    if [ "\$N_BINS" -gt 0 ]; then
        export GTDBTK_DATA_PATH=${params.db_root}/gtdbtk_r232/release232

        gtdbtk classify_wf \
            --genome_dir ${bins_dir} \
            --out_dir gtdbtk_results \
            --extension fa \
            --cpus ${task.cpus} \
            --pplacer_cpus 1 \
            --scratch_dir \${TMPDIR:-/tmp}/gtdbtk_${sample}

        # Summary: extract classifications
        cat gtdbtk_results/classify/gtdbtk.*.summary.tsv > ${sample}_taxonomy.tsv 2>/dev/null || true

        rm -rf \${TMPDIR:-/tmp}/gtdbtk_${sample}
    else
        echo "No bins for ${sample}, skipping GTDB-Tk"
        mkdir -p gtdbtk_results
        echo -e "user_genome\tclassification" > ${sample}_taxonomy.tsv
    fi
    """
}


// ======================================================================
// STEP 10b: SOURMASH_CLASSIFY — ultra-light GTDB-Tk alternative
// Used when params.taxonomy_tool == 'sourmash'. K-mer sketching against
// the prepared sourmash GTDB database (~3.7 GB). RAM ≤8 GB, runtime
// seconds-to-minutes per MAG.
// ======================================================================

process SOURMASH_CLASSIFY {
    tag "${sample}"
    publishDir "${outdir}/24_taxonomy_sourmash/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bins_dir)

    output:
    path("${sample}_taxonomy.tsv"), emit: summary
    path("gather_results/"),        emit: gather, optional: true

    script:
    def db_path = "${params.db_root}/sourmash/gtdb-rs226-reps-k31.zip"
    def lineages = "${params.db_root}/sourmash/gtdb-rs226.lineages.csv"
    """
    mkdir -p gather_results sketches

    N_BINS=\$(ls ${bins_dir}/*.fa 2>/dev/null | wc -l)

    if [ "\$N_BINS" -gt 0 ] && [ -f "${db_path}" ]; then
        # 1. Sketch cada MAG (k=31, scaled=1000)
        for fa in ${bins_dir}/*.fa; do
            bin_name=\$(basename "\$fa" .fa)
            sourmash sketch dna -p k=31,scaled=1000 "\$fa" \
                -o "sketches/\${bin_name}.sig" \
                --name "\${bin_name}" 2>/dev/null || true
        done

        # 2. Gather contra GTDB
        for sig in sketches/*.sig; do
            bin_name=\$(basename "\$sig" .sig)
            sourmash gather "\$sig" "${db_path}" \
                -o "gather_results/\${bin_name}.gather.csv" \
                --threshold-bp 50000 \
                -k 31 2>/dev/null || true
        done

        # 3. Concatenate gather results into a single CSV.
        # Bug guard: `head -1 gather_results/*.csv` with multiple files prints
        # "==> file <==" banners that would corrupt the header. Pick one file
        # explicitly to extract the header instead.
        first_gather=\$(ls gather_results/*.gather.csv 2>/dev/null | head -1)
        if [ -n "\$first_gather" ]; then
            head -1 "\$first_gather" > gather_all.csv
            for f in gather_results/*.gather.csv; do
                tail -n +2 "\$f" >> gather_all.csv
            done
        else
            : > gather_all.csv
        fi

        if [ -f "${lineages}" ] && [ -s gather_all.csv ]; then
            sourmash tax annotate -g gather_all.csv \
                -t "${lineages}" \
                -o tax_annotated 2>/dev/null || true

            # Also generate per-MAG classification
            sourmash tax genome -g gather_all.csv \
                -t "${lineages}" \
                --output-format csv_summary \
                -o ${sample}_sourmash_tax.csv 2>/dev/null || true
        fi

        # 4. Emit a TSV compatible with the GTDB-Tk format
        echo -e "user_genome\\tclassification\\tani\\tcontainment" > ${sample}_taxonomy.tsv
        if [ -f "${sample}_sourmash_tax.csv" ]; then
            # sourmash tax genome output: query_name,status,rank,fraction,lineage,...
            tail -n +2 ${sample}_sourmash_tax.csv | while IFS=',' read -r qname status rank fraction lineage rest; do
                echo -e "\${qname}\\t\${lineage}\\t\\t\${fraction}" >> ${sample}_taxonomy.tsv
            done
        fi
    else
        echo "No bins or sourmash DB not found"
        echo -e "user_genome\\tclassification\\tani\\tcontainment" > ${sample}_taxonomy.tsv
    fi
    """
}


// ======================================================================
// STEP 10c: SKANI_CLASSIFY — fast ANI-based MAG classification (default)
// Used when params.taxonomy_tool == 'skani'. Real ANI against the
// pre-built skani GTDB r226 sketch DB (~36 GB, v0.3 single-file layout).
// RAM ≤16-24 GB, runtime seconds per MAG. Closest GTDB-Tk replacement
// for resource-constrained environments.
// ======================================================================

process SKANI_CLASSIFY {
    tag "${sample}"
    publishDir "${outdir}/24_taxonomy_skani/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bins_dir)

    output:
    path("${sample}_taxonomy.tsv"),  emit: summary
    path("${sample}_skani_full.tsv"), emit: full, optional: true

    script:
    def db_path = "${params.db_root}/skani/skani_gtdb_r226-v0.3"
    """
    set -euo pipefail
    shopt -s nullglob

    # Always emit a header even when there are no bins, so downstream
    # steps that expect the file (and the parsing code in
    # AMR_PATHOGEN_INTEGRATION) never trip on a missing file.
    echo -e "user_genome\\tclassification\\tani\\taf" > ${sample}_taxonomy.tsv

    bins=( ${bins_dir}/*.fa )
    if [ \${#bins[@]} -eq 0 ] || [ ! -d "${db_path}" ]; then
        echo "[skani] No bins or skani DB not found at ${db_path}, emitting empty taxonomy"
        exit 0
    fi

    # `skani search -d <db_dir>` runs each query against the sketched DB.
    # The default output table includes ANI, Align_fraction_query and
    # Ref_name (which embeds the GTDB lineage in the sketch metadata).
    skani search \\
        -d ${db_path} \\
        -q "\${bins[@]}" \\
        -t ${task.cpus} \\
        -o ${sample}_skani_full.tsv

    # Reduce to one best-ANI hit per query and write a GTDB-Tk-shaped TSV.
    python3 - <<'PYEOF'
    import csv, os, sys
    inp  = "${sample}_skani_full.tsv"
    outp = "${sample}_taxonomy.tsv"
    best = {}
    with open(inp) as fh:
        reader = csv.DictReader(fh, delimiter='\\t')
        for row in reader:
            q_path = row.get('Query_file') or row.get('Query_name') or ''
            q = os.path.basename(q_path)
            if q.endswith('.fa'):
                q = q[:-3]
            try:
                ani = float(row.get('ANI', '') or 0)
            except ValueError:
                ani = 0.0
            try:
                af = float(row.get('Align_fraction_query', '') or 0)
            except ValueError:
                af = 0.0
            ref = row.get('Ref_name', '') or row.get('Ref_file', '')
            if not q:
                continue
            if q not in best or ani > best[q][1]:
                best[q] = (ref, ani, af)
    with open(outp, 'a') as fh:
        for q, (ref, ani, af) in best.items():
            fh.write(f"{q}\\t{ref}\\t{ani:.2f}\\t{af:.3f}\\n")
    PYEOF
    """
}


// ======================================================================
// STEP 11: AMRFINDERPLUS — AMR + virulence detection in MAGs
// ======================================================================

process AMRFINDERPLUS_MAGS {
    tag "${sample}"
    publishDir "${outdir}/25_amr_mags/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bins_dir)

    output:
    path("${sample}_amr_combined.tsv"), emit: report
    path("per_bin/"),                   emit: per_bin

    script:
    """
    mkdir -p per_bin

    N_BINS=\$(ls ${bins_dir}/*.fa 2>/dev/null | wc -l)

    if [ "\$N_BINS" -gt 0 ]; then
        for bin_fa in ${bins_dir}/*.fa; do
            bin_name=\$(basename "\$bin_fa" .fa)

            # Capture stderr to per-bin log so failures are visible (the
            # earlier `2>/dev/null` was swallowing the "BLAST DB not
            # found" message that hid a wrong DB path for hours).
            amrfinder \
                --nucleotide "\$bin_fa" \
                --database ${params.amrfinder_db} \
                --threads ${task.cpus} \
                --plus \
                --name "\$bin_name" \
                --output "per_bin/\${bin_name}_amr.tsv" \
                2>"per_bin/\${bin_name}.log" || true
        done

        # Concatenate per-bin TSVs. Bug guard: `head -1 *.tsv` over a
        # multi-file glob prints `==> file <==` banners; explicitly read
        # the header from the first file instead.
        first_tsv=\$(ls per_bin/*_amr.tsv 2>/dev/null | head -1)
        if [ -n "\$first_tsv" ]; then
            head -1 "\$first_tsv" > ${sample}_amr_combined.tsv
            for f in per_bin/*_amr.tsv; do
                tail -n +2 "\$f" >> ${sample}_amr_combined.tsv
            done
        else
            : > ${sample}_amr_combined.tsv
        fi
    else
        echo "No bins for ${sample}"
        echo -e "Name\tProtein identifier\tContig id\tStart\tStop\tStrand\tGene symbol\tSequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference sequence\t% Identity to reference sequence\tAlignment length\tAccession of closest sequence\tName of closest sequence\tHMM id\tHMM description" > ${sample}_amr_combined.tsv
        mkdir -p per_bin
    fi
    """
}


// ======================================================================
// STEP 12: GENOMAD — plasmid + virus detection on assemblies
// ======================================================================

process GENOMAD {
    tag "${sample}"
    publishDir "${outdir}/26_genomad/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly)

    output:
    path("genomad_output/"),                                       emit: results
    path("${sample}_plasmid_summary.tsv"),                         emit: plasmid_summary, optional: true
    path("${sample}_virus_summary.tsv"),                           emit: virus_summary, optional: true
    path("${sample}_plasmid_genes.tsv"),                           emit: plasmid_genes, optional: true
    tuple val(sample), path("${sample}_plasmid.fna"),              emit: plasmid_fna, optional: true

    script:
    """
    genomad end-to-end \
        ${assembly} \
        genomad_output \
        ${params.db_root}/genomad/genomad_db \
        --threads ${task.cpus}

    # Extract plasmid summary
    if [ -f genomad_output/*_summary/*_plasmid_summary.tsv ]; then
        cp genomad_output/*_summary/*_plasmid_summary.tsv ${sample}_plasmid_summary.tsv
    fi

    # Extract virus summary
    if [ -f genomad_output/*_summary/*_virus_summary.tsv ]; then
        cp genomad_output/*_summary/*_virus_summary.tsv ${sample}_virus_summary.tsv
    fi

    # Genes on plasmids (used to cross-reference with AMR)
    if [ -f genomad_output/*_summary/*_plasmid_genes.tsv ]; then
        cp genomad_output/*_summary/*_plasmid_genes.tsv ${sample}_plasmid_genes.tsv
    fi

    # Plasmid contigs as FASTA, for an independent AMRFinderPlus pass.
    if [ -f genomad_output/*_summary/*_plasmid.fna ]; then
        cp genomad_output/*_summary/*_plasmid.fna ${sample}_plasmid.fna
    fi
    """
}


// ======================================================================
// STEP 12b: AMR_ON_PLASMIDS — AMRFinderPlus on geNomad plasmid contigs
// Independent of AMR-in-MAGs; surfaces resistance markers carried by
// contigs that may transfer horizontally regardless of MAG assignment.
// ======================================================================

process AMR_ON_PLASMIDS {
    tag "${sample}"
    publishDir "${outdir}/26_genomad/${sample}", mode: 'copy', pattern: "${sample}_plasmid_amr.tsv"

    input:
    tuple val(sample), path(plasmid_fna)

    output:
    path("${sample}_plasmid_amr.tsv"), emit: tsv, optional: true

    script:
    """
    # Empty FASTA => emit a header-only TSV so downstream steps stay tidy.
    if [ ! -s "${plasmid_fna}" ]; then
        printf 'Sample\tContig\tGene\tName\tClass\tSubclass\tMethod\tIdentity\tCoverage\tAccession\tClosest_ref\\n' \
            > "${sample}_plasmid_amr.tsv"
        exit 0
    fi

    amrfinder \\
        --nucleotide "${plasmid_fna}" \\
        --database ${params.amrfinder_db} \\
        --threads ${task.cpus} \\
        --plus \\
        --name "${sample}" \\
        --output "${sample}_plasmid_amr_raw.tsv" \\
        2>"${sample}_plasmid_amr.log" || true

    # Re-shape AMRFinderPlus output into the columns the report expects.
    python3 - "${sample}" "${sample}_plasmid_amr_raw.tsv" \\
        > "${sample}_plasmid_amr.tsv" <<'PY'
import sys, csv
sample, src = sys.argv[1], sys.argv[2]
out_cols = ['Sample','Contig','Gene','Name','Class','Subclass','Method',
            'Identity','Coverage','Accession','Closest_ref']
print('\\t'.join(out_cols))
try:
    with open(src) as fh:
        reader = csv.DictReader(fh, delimiter='\\t')
        for r in reader:
            print('\\t'.join([
                sample,
                r.get('Contig id', ''),
                r.get('Element symbol', r.get('Gene symbol', '')),
                r.get('Element name', r.get('Sequence name', '')),
                r.get('Class', ''),
                r.get('Subclass', ''),
                r.get('Method', ''),
                r.get('% Identity to reference',
                      r.get('% Identity to reference sequence', '')),
                r.get('% Coverage of reference',
                      r.get('% Coverage of reference sequence', '')),
                r.get('Closest reference accession', ''),
                r.get('Closest reference name', ''),
            ]))
except FileNotFoundError:
    pass
PY
    """
}


// ======================================================================
// STEP 13: ABRICATE — virulence (VFDB) + AMR (CARD) + PlasmidFinder
// ======================================================================

process ABRICATE {
    tag "${sample}"
    publishDir "${outdir}/28_abricate/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bins_dir)

    output:
    path("*.vfdb.tsv"),          emit: vfdb, optional: true
    path("*.card.tsv"),          emit: card, optional: true
    path("*.plasmidfinder.tsv"), emit: plasmidfinder, optional: true

    script:
    """
    N_BINS=\$(ls ${bins_dir}/*.fa 2>/dev/null | wc -l)
    if [ "\$N_BINS" -gt 0 ]; then
        for bin_fa in ${bins_dir}/*.fa; do
            bin_name=\$(basename "\$bin_fa" .fa)
            for db in vfdb card plasmidfinder; do
                abricate --db "\$db" --threads ${task.cpus} "\$bin_fa" \
                    > "\${bin_name}.\${db}.tsv" 2>/dev/null || true
            done
        done
    fi
    """
}


// ======================================================================
// STEP 14: MOBSUITE — plasmid typing
// ======================================================================

process MOBSUITE {
    tag "${sample}"
    publishDir "${outdir}/29_mobsuite/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(assembly)

    output:
    path("mobtyper_results.txt"), emit: results, optional: true
    path("*"),                    emit: all_files

    script:
    """
    mob_recon --infile ${assembly} --outdir . \
        --num_threads ${task.cpus} --force 2>/dev/null || true
    touch mobtyper_results.txt
    """
}


// ======================================================================
// STEP 15: INTEGRONFINDER — integron detection per MAG
// ======================================================================

process INTEGRONFINDER {
    tag "${sample}/${bin_name}"
    publishDir "${outdir}/30_integronfinder/${sample}/${bin_name}", mode: 'copy'

    input:
    tuple val(sample), val(bin_name), path(bin_fa)

    output:
    path("Results_*/*.integrons"), emit: integrons, optional: true
    path("Results_*/*.summary"),   emit: summary, optional: true

    script:
    """
    integron_finder --local-max --func-annot \
        --cpu ${task.cpus} --outdir . \
        ${bin_fa} 2>/dev/null || true
    """
}


// ======================================================================
// STEP 16: BAKTA — functional annotation per MAG
// ======================================================================

process BAKTA {
    tag "${sample}/${bin_name}"
    publishDir "${outdir}/31_bakta/${sample}/${bin_name}", mode: 'copy'

    input:
    tuple val(sample), val(bin_name), path(bin_fa)

    output:
    path("${bin_name}.gff3"),  emit: gff, optional: true
    path("${bin_name}.tsv"),   emit: tsv, optional: true
    path("${bin_name}.txt"),   emit: summary, optional: true
    path("${bin_name}.*"),     emit: all_files

    script:
    def locus = sample.replaceAll('[^a-zA-Z0-9_-]', '_').take(24)
    """
    export TMPDIR=\${TMPDIR:-/tmp}
    bakta --db ${params.bakta_db} \
        --output . --prefix ${bin_name} \
        --threads ${task.cpus} \
        --locus-tag ${locus} \
        --skip-plot --force \
        ${bin_fa} 2>${bin_name}.bakta.log || true
    """
}


// ======================================================================
// STEP 17a: AMR_PATHOGEN_INTEGRATION — combine AMR + taxonomy + plasmids
// ======================================================================

process AMR_PATHOGEN_INTEGRATION {
    tag "${sample}"
    publishDir "${outdir}/27_amr_pathogen_integration", mode: 'copy'

    input:
    tuple val(sample), path(taxonomy), path(amr), path(checkm2_qc), path(plasmids)

    output:
    path("${sample}_amr_pathogen.tsv"),     emit: report
    path("${sample}_summary.txt"),          emit: summary

    script:
    """
    python3 ${params.scripts_dir}/integrate_amr_pathogen.py \
        --sample ${sample} \
        --taxonomy ${taxonomy} \
        --amr ${amr} \
        --checkm2 ${checkm2_qc} \
        --plasmids ${plasmids} \
        --output ${sample}_amr_pathogen.tsv \
        --summary ${sample}_summary.txt
    """
}


// ======================================================================
// STEP 17b: MAG_REPORT — final HTML + Excel report
// ======================================================================

process MAG_REPORT {
    publishDir "${outdir}/32_reports", mode: 'copy'

    // The MAG report depends on every upstream MAG-phase artefact, but
    // it reads them from `${outdir}` directly — the input channel is
    // only a synchronisation barrier. We pass the trigger paths as
    // `val` (string list) instead of `path` so Nextflow does not stage
    // them into this workdir; staging would collide because the per-bin
    // outputs of Bakta/IntegronFinder share basenames across samples
    // (every sample has its own SemiBin_0.fna, SemiBin_0.gbff …) and
    // Nextflow refuses to stage two files with the same name.
    // Cache invalidation still works: the val list participates in the
    // task hash, so any change in upstream paths re-runs the report.
    input:
    val(triggers)

    output:
    path("*.html"), emit: html
    path("*.xlsx"), emit: excel, optional: true

    script:
    """
    python3 ${params.scripts_dir}/generate_mag_report_v2.py \
        --results-dir ${outdir} \
        --run-name ${run_name} \
        --input-dir ${params.input} \
        --branding ${projectDir}/branding/ \
        --output ${run_name}_epitaxmag_report.html || true

    python3 ${params.scripts_dir}/generate_excel_report.py \
        --results-dir ${outdir} \
        --run-name ${run_name} \
        --output ${run_name}_full_report.xlsx || true
    """
}


// ======================================================================
// WORKFLOW MAG
// ======================================================================

workflow MAG {

    take:
    ch_filtered    // tuple(sample, fastq) — reads filtered by TAX

    main:

    // -- Assembly (skippable with --assemblies_dir)
    if (params.assemblies_dir) {
        log.info "  MAG: reusing existing assemblies from ${params.assemblies_dir}"
        ch_assemblies = Channel
            .fromPath("${params.assemblies_dir}/*/*.assembly.fasta", checkIfExists: true)
            .map { f -> tuple(f.parent.name, f) }
    } else {
        METAFLYE(ch_filtered)
        ch_assemblies = METAFLYE.out.assembly
    }

    // -- Polishing
    ch_medaka_input = ch_assemblies
        .join(ch_filtered)
        .map { sample, assembly, fastq -> tuple(sample, assembly, fastq) }
    MEDAKA(ch_medaka_input)

    // -- Assembly QC
    ch_quast_assemblies = MEDAKA.out.polished.map { s, a -> a }.collect()
    ch_quast_labels = MEDAKA.out.polished.map { s, a -> s }.collect()
    QUAST(ch_quast_assemblies, ch_quast_labels)

    // -- Coverage mapping
    ch_map_input = MEDAKA.out.polished
        .join(ch_filtered)
        .map { sample, assembly, fastq -> tuple(sample, assembly, fastq) }
    MAP_READS(ch_map_input)

    // -- Binning (3 binners in parallel)
    ch_bin_input = MAP_READS.out.assembly
        .join(MAP_READS.out.bam)
        .map { sample, assembly, bam, bai -> tuple(sample, assembly, bam, bai) }
    METABAT2(ch_bin_input)
    SEMIBIN2(ch_bin_input)
    ch_maxbin_input = MEDAKA.out.polished
        .join(ch_filtered)
        .map { sample, assembly, fastq -> tuple(sample, assembly, fastq) }
    MAXBIN2(ch_maxbin_input)

    // -- Refinement
    ch_dastool_input = MAP_READS.out.assembly
        .join(METABAT2.out.bins)
        .join(MAXBIN2.out.bins)
        .join(SEMIBIN2.out.bins)
        .map { sample, assembly, metabat, maxbin, semibin ->
            tuple(sample, assembly, metabat, maxbin, semibin) }
    DAS_TOOL(ch_dastool_input)

    // -- MAG QC
    CHECKM2(DAS_TOOL.out.refined_bins)

    // -- Taxonomy (selected by params.taxonomy_tool).
    //   skani    : default. Real ANI, ~16-24 GB RAM, ~36 GB DB (GTDB r226).
    //   sourmash : ultra-light fallback, ~5 GB RAM, ~3.7 GB DB.
    //   gtdbtk   : gold-standard reference, ≥120 GB RAM, ~60 GB DB.
    switch (params.taxonomy_tool) {
        case 'skani':
            SKANI_CLASSIFY(DAS_TOOL.out.refined_bins)
            break
        case 'sourmash':
            log.info "  Using Sourmash for MAG taxonomy (~5 GB RAM, k-mer based)."
            SOURMASH_CLASSIFY(DAS_TOOL.out.refined_bins)
            break
        case 'gtdbtk':
            log.info "  Using GTDB-Tk for MAG taxonomy (≥120 GB RAM, gold standard)."
            GTDBTK(DAS_TOOL.out.refined_bins)
            break
        default:
            error "Unknown taxonomy_tool '${params.taxonomy_tool}'. Use 'skani', 'sourmash' or 'gtdbtk'."
    }

    // -- AMR on MAGs
    AMRFINDERPLUS_MAGS(DAS_TOOL.out.refined_bins)

    // -- Plasmids and viruses (on the full assembly)
    GENOMAD(MEDAKA.out.polished)

    // -- AMRFinderPlus on the plasmid contigs predicted by geNomad.
    // Independent of AMR-in-MAGs; useful for AMR carried on contigs that
    // never made it into a high-quality bin.
    AMR_ON_PLASMIDS(GENOMAD.out.plasmid_fna)

    // -- Virulence (VFDB) + AMR (CARD) + PlasmidFinder
    ABRICATE(DAS_TOOL.out.refined_bins)

    // -- Plasmid typing
    MOBSUITE(MEDAKA.out.polished)

    // -- Integrons (per MAG)
    ch_integronfinder_input = DAS_TOOL.out.refined_bins
        .flatMap { sample, bins_dir ->
            def fa_files = file("${bins_dir}/*.fa")
            fa_files.collect { fa ->
                tuple(sample, fa.baseName, fa)
            }
        }
    INTEGRONFINDER(ch_integronfinder_input)

    // -- Bakta annotation (per MAG)
    ch_bakta_input = DAS_TOOL.out.refined_bins
        .flatMap { sample, bins_dir ->
            def fa_files = file("${bins_dir}/*.fa")
            fa_files.collect { fa ->
                tuple(sample, fa.baseName, fa)
            }
        }
    BAKTA(ch_bakta_input)

    // -- AMR + pathogen + plasmid integration (per sample).
    // The upstream processes emit bare paths whose filenames carry the
    // sample id (e.g. "<sample>_taxonomy.tsv"). We rebuild
    // `tuple(sample, file)` channels here via map+filename so we can
    // join across taxonomy, AMR, CheckM2 and plasmid summary without
    // changing the upstream outputs (which would invalidate caches).
    // Uses the output of whichever MAG taxonomy tool was selected.
    ch_tax_file = (params.taxonomy_tool == 'skani'    ? SKANI_CLASSIFY.out.summary    :
                   params.taxonomy_tool == 'sourmash' ? SOURMASH_CLASSIFY.out.summary :
                                                        GTDBTK.out.summary)
    ch_taxonomy = ch_tax_file
        .map { f -> tuple(f.baseName.replaceAll(/_taxonomy$/, ''), f) }
    ch_amr_pathogen_join = AMRFINDERPLUS_MAGS.out.report
        .map { f -> tuple(f.baseName.replaceAll(/_amr_combined$/, ''), f) }
    // CHECKM2 already emits a (sample, file) tuple, so it can be joined
    // directly without filename parsing.
    ch_checkm2_pathogen = CHECKM2.out.report
    ch_plasmids_pathogen = GENOMAD.out.plasmid_summary
        .map { f -> tuple(f.baseName.replaceAll(/_plasmid_summary$/, ''), f) }

    ch_pathogen_input = ch_taxonomy
        .join(ch_amr_pathogen_join)
        .join(ch_checkm2_pathogen)
        .join(ch_plasmids_pathogen)
    AMR_PATHOGEN_INTEGRATION(ch_pathogen_input)

    // -- Final HTML/Excel report. Use the actual artefact paths from the
    // upstream per-sample / per-bin steps as inputs so Nextflow's cache
    // tracks them properly: when any of these regenerate, MAG_REPORT
    // re-runs instead of serving a stale cached HTML.
    // Mix the artefact paths from every upstream producer; the mix
    // can include lists (per-bin outputs) so we flatten and stringify
    // before .collect() — `val` items must be Strings/serialisable for
    // the task hash to be stable across runs.
    ch_report_inputs = AMR_PATHOGEN_INTEGRATION.out.report
        .mix(
            BAKTA.out.all_files,
            INTEGRONFINDER.out.summary.ifEmpty([]),
            ABRICATE.out.vfdb.ifEmpty([]),
            MOBSUITE.out.results.ifEmpty([]),
            AMR_ON_PLASMIDS.out.tsv.ifEmpty([])
        )
        .flatten()
        .map { it.toString() }
        .collect()
    MAG_REPORT(ch_report_inputs)

    emit:
    assemblies    = MEDAKA.out.polished
    refined_bins  = DAS_TOOL.out.refined_bins
    checkm2       = CHECKM2.out.report
    taxonomy      = ch_tax_file
    amr_mags      = AMRFINDERPLUS_MAGS.out.report
    plasmids      = GENOMAD.out.plasmid_summary
    report        = MAG_REPORT.out.html
}
