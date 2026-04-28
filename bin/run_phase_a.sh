#!/bin/bash
# EpiTaxMAG Fase A — Modulos adicionales sobre resultados existentes
# Sin tocar Nextflow. Corre directamente con Singularity.
#
# Uso:
#   bash bin/run_phase_a.sh results/260226_EPIM232
#
# Requisitos:
#   - Singularity instalado
#   - Resultados de MAG completados (pasos 14-25)
#   - Containers se descargan automaticamente al cache

set -euo pipefail

RESULTS="${1:?Uso: bash bin/run_phase_a.sh results/<run_name>}"
RUN_NAME=$(basename "$RESULTS")
THREADS="${2:-16}"

CACHE="${NXF_SINGULARITY_CACHEDIR:-/GRU0/DATABASES/singularity_cache}"
BIND="${SINGULARITY_BIND:---bind /GRU0 --bind /home/asanzc/jobs}"
DB_ROOT="${DB_ROOT:-/home/asanzc/jobs/nanometa/resources}"

# Containers
SIF_ABRICATE="$CACHE/depot.galaxyproject.org-singularity-abricate%3A1.4.0--h05cac1d_0.img"
SIF_MOBSUITE="$CACHE/depot.galaxyproject.org-singularity-mob_suite%3A3.1.9--pyhdfd78af_1.img"
SIF_BAKTA="$CACHE/depot.galaxyproject.org-singularity-bakta%3A1.12.0--pyhdfd78af_0.img"
SIF_INTFINDER="$CACHE/depot.galaxyproject.org-singularity-integron_finder%3A2.0rc6--py_0.img"
SIF_AMRFINDER="$CACHE/depot.galaxyproject.org-singularity-ncbi-amrfinderplus%3A4.2.7--hf69ffd2_0.img"

pull_if_missing() {
    local sif="$1" url="$2"
    if [ ! -f "$sif" ]; then
        echo "  Descargando $(basename "$sif")..."
        singularity pull --name "$(basename "$sif")" "$url"
        mv "$(basename "$sif")" "$sif"
    fi
}

srun() {
    singularity exec --no-home $BIND "$@"
}

echo ""
echo "  EpiTaxMAG Fase A — Modulos adicionales"
echo "  Run: $RUN_NAME"
echo "  Threads: $THREADS"
echo "  Cache: $CACHE"
echo ""

# Pull containers si faltan
echo "[1/6] Verificando containers..."
pull_if_missing "$SIF_ABRICATE" "https://depot.galaxyproject.org/singularity/abricate%3A1.4.0--h05cac1d_0"
pull_if_missing "$SIF_MOBSUITE" "https://depot.galaxyproject.org/singularity/mob_suite%3A3.1.9--pyhdfd78af_1"
pull_if_missing "$SIF_BAKTA" "https://depot.galaxyproject.org/singularity/bakta%3A1.12.0--pyhdfd78af_0"
pull_if_missing "$SIF_INTFINDER" "https://depot.galaxyproject.org/singularity/integron_finder%3A2.0rc6--py_0"
echo "  Containers OK"

# ══════════════════════════════════════════════════════════════
# MODULO 1: ABRicate — Virulencia (VFDB) + PlasmidFinder
# ══════════════════════════════════════════════════════════════
ABRICATE_DIR="$RESULTS/30_abricate"
echo ""
echo "[2/6] ABRicate (VFDB + CARD + PlasmidFinder)..."

if [ -d "$ABRICATE_DIR" ] && [ "$(find "$ABRICATE_DIR" -name "*.vfdb.tsv" 2>/dev/null | wc -l)" -gt 0 ]; then
    echo "  Ya existe, saltando"
else
    mkdir -p "$ABRICATE_DIR"
    for sample_dir in "$RESULTS"/21_binning_dastool/*/; do
        sample=$(basename "$sample_dir")
        bins_dir="$sample_dir/${sample}_DASTool_bins"
        [ ! -d "$bins_dir" ] && continue
        n_bins=$(ls "$bins_dir"/*.fa 2>/dev/null | wc -l)
        [ "$n_bins" -eq 0 ] && continue

        mkdir -p "$ABRICATE_DIR/$sample"

        for bin_fa in "$bins_dir"/*.fa; do
            bin_name=$(basename "$bin_fa" .fa)
            for db in vfdb card plasmidfinder; do
                out="$ABRICATE_DIR/$sample/${bin_name}.${db}.tsv"
                [ -f "$out" ] && continue
                srun "$SIF_ABRICATE" abricate --db "$db" --threads "$THREADS" "$bin_fa" > "$out" 2>/dev/null || true
            done
        done

        # Resumen por muestra
        srun "$SIF_ABRICATE" abricate --summary "$ABRICATE_DIR/$sample"/*.vfdb.tsv \
            > "$ABRICATE_DIR/$sample/${sample}_vfdb_summary.tsv" 2>/dev/null || true
        srun "$SIF_ABRICATE" abricate --summary "$ABRICATE_DIR/$sample"/*.card.tsv \
            > "$ABRICATE_DIR/$sample/${sample}_card_summary.tsv" 2>/dev/null || true

        echo "  $sample: $n_bins MAGs procesados"
    done
fi

# ══════════════════════════════════════════════════════════════
# MODULO 2: MOB-suite — Caracterizacion de plasmidos
# ══════════════════════════════════════════════════════════════
MOBSUITE_DIR="$RESULTS/31_mobsuite"
echo ""
echo "[3/6] MOB-suite (tipificacion de plasmidos)..."

if [ -d "$MOBSUITE_DIR" ] && [ "$(find "$MOBSUITE_DIR" -name "mobtyper_results.txt" 2>/dev/null | wc -l)" -gt 0 ]; then
    echo "  Ya existe, saltando"
else
    mkdir -p "$MOBSUITE_DIR"
    for asm in "$RESULTS"/15_polish_medaka/*.polished.fasta; do
        [ ! -f "$asm" ] && continue
        sample=$(basename "$asm" .polished.fasta)
        out_dir="$MOBSUITE_DIR/$sample"
        [ -d "$out_dir" ] && [ -f "$out_dir/mobtyper_results.txt" ] && continue

        mkdir -p "$out_dir"
        echo "  $sample: corriendo mob_recon..."
        srun "$SIF_MOBSUITE" mob_recon \
            --infile "$asm" \
            --outdir "$out_dir" \
            --num_threads "$THREADS" \
            --force 2>/dev/null || true
    done
fi

# ══════════════════════════════════════════════════════════════
# MODULO 3: IntegronFinder — Integrones
# ══════════════════════════════════════════════════════════════
INTFINDER_DIR="$RESULTS/32_integronfinder"
echo ""
echo "[4/6] IntegronFinder (integrones de clase 1/2/3)..."

if [ -d "$INTFINDER_DIR" ] && [ "$(find "$INTFINDER_DIR" -name "*.integrons" 2>/dev/null | wc -l)" -gt 0 ]; then
    echo "  Ya existe, saltando"
else
    mkdir -p "$INTFINDER_DIR"
    for asm in "$RESULTS"/15_polish_medaka/*.polished.fasta; do
        [ ! -f "$asm" ] && continue
        sample=$(basename "$asm" .polished.fasta)
        out_dir="$INTFINDER_DIR/$sample"
        [ -d "$out_dir" ] && continue

        echo "  $sample: buscando integrones..."
        # IntegronFinder needs writable tmp
        export TMPDIR="$RESULTS/../work/tmp"
        mkdir -p "$TMPDIR"
        srun "$SIF_INTFINDER" integron_finder \
            --local-max --func-annot \
            --cpu "$THREADS" \
            --outdir "$out_dir" \
            "$asm" 2>/dev/null || true
    done
fi

# ══════════════════════════════════════════════════════════════
# MODULO 4: AMRFinderPlus con --organism (mutaciones puntuales)
# ══════════════════════════════════════════════════════════════
AMR_ORG_DIR="$RESULTS/33_amr_organism"
echo ""
echo "[5/6] AMRFinderPlus --organism (mutaciones puntuales)..."

if [ -d "$AMR_ORG_DIR" ] && [ "$(find "$AMR_ORG_DIR" -name "*_amr_org.tsv" 2>/dev/null | wc -l)" -gt 0 ]; then
    echo "  Ya existe, saltando"
else
    mkdir -p "$AMR_ORG_DIR"

    # Map GTDB genus to AMRFinderPlus organism names
    declare -A GENUS_TO_ORG=(
        ["Acinetobacter"]="Acinetobacter_baumannii"
        ["Escherichia"]="Escherichia"
        ["Klebsiella"]="Klebsiella_pneumoniae"
        ["Pseudomonas"]="Pseudomonas_aeruginosa"
        ["Salmonella"]="Salmonella"
        ["Staphylococcus"]="Staphylococcus_aureus"
        ["Enterococcus"]="Enterococcus_faecalis"
        ["Neisseria"]="Neisseria_meningitidis"
        ["Campylobacter"]="Campylobacter"
        ["Clostridioides"]="Clostridioides_difficile"
        ["Burkholderia"]="Burkholderia_cepacia"
        ["Vibrio"]="Vibrio_cholerae"
    )

    for sample_dir in "$RESULTS"/21_binning_dastool/*/; do
        sample=$(basename "$sample_dir")
        bins_dir="$sample_dir/${sample}_DASTool_bins"
        tax_file="$RESULTS/23_taxonomy_gtdbtk/$sample/${sample}_taxonomy.tsv"

        [ ! -d "$bins_dir" ] && continue
        [ ! -f "$tax_file" ] && continue

        mkdir -p "$AMR_ORG_DIR/$sample"

        # Parse taxonomy to get genus per bin
        while IFS=$'\t' read -r bin_name taxonomy rest; do
            [ "$bin_name" = "user_genome" ] && continue
            bin_fa="$bins_dir/${bin_name}.fa"
            [ ! -f "$bin_fa" ] && continue

            # Extract genus from GTDB taxonomy
            genus=$(echo "$taxonomy" | grep -oP 'g__\K[^;]+' || echo "")
            organism="${GENUS_TO_ORG[$genus]:-}"

            out="$AMR_ORG_DIR/$sample/${bin_name}_amr_org.tsv"
            [ -f "$out" ] && continue

            if [ -n "$organism" ]; then
                srun "$SIF_AMRFINDER" amrfinder \
                    --nucleotide "$bin_fa" \
                    --organism "$organism" \
                    --database "$DB_ROOT/amrfinderplus/latest" \
                    --threads "$THREADS" \
                    --plus \
                    --name "$bin_name" \
                    --output "$out" 2>/dev/null || true
                echo "  $bin_name ($genus -> $organism)"
            else
                # No organism mapping, run without --organism
                srun "$SIF_AMRFINDER" amrfinder \
                    --nucleotide "$bin_fa" \
                    --database "$DB_ROOT/amrfinderplus/latest" \
                    --threads "$THREADS" \
                    --plus \
                    --name "$bin_name" \
                    --output "$out" 2>/dev/null || true
            fi
        done < "$tax_file"
    done
fi

# ══════════════════════════════════════════════════════════════
# MODULO 5: Bakta — Anotacion funcional de MAGs
# ══════════════════════════════════════════════════════════════
BAKTA_DIR="$RESULTS/34_bakta"
echo ""
echo "[6/6] Bakta (anotacion funcional de MAGs)..."

if [ -d "$BAKTA_DIR" ] && [ "$(find "$BAKTA_DIR" -name "*.gff3" 2>/dev/null | wc -l)" -gt 0 ]; then
    echo "  Ya existe, saltando"
else
    mkdir -p "$BAKTA_DIR"
    BAKTA_DB="$DB_ROOT/bakta/db"

    if [ ! -d "$BAKTA_DB" ]; then
        echo "  AVISO: Bakta DB no encontrada en $BAKTA_DB. Saltando Bakta."
    else
        for sample_dir in "$RESULTS"/21_binning_dastool/*/; do
            sample=$(basename "$sample_dir")
            bins_dir="$sample_dir/${sample}_DASTool_bins"
            [ ! -d "$bins_dir" ] && continue
            n_bins=$(ls "$bins_dir"/*.fa 2>/dev/null | wc -l)
            [ "$n_bins" -eq 0 ] && continue

            mkdir -p "$BAKTA_DIR/$sample"

            for bin_fa in "$bins_dir"/*.fa; do
                bin_name=$(basename "$bin_fa" .fa)
                out_dir="$BAKTA_DIR/$sample/$bin_name"
                [ -d "$out_dir" ] && [ -f "$out_dir/${bin_name}.gff3" ] && continue

                echo "  $bin_name..."
                srun "$SIF_BAKTA" bakta \
                    --db "$BAKTA_DB" \
                    --output "$out_dir" \
                    --prefix "$bin_name" \
                    --threads "$THREADS" \
                    --locus-tag "$bin_name" \
                    --skip-plot \
                    "$bin_fa" 2>/dev/null || true
            done
            echo "  $sample: $n_bins MAGs anotados"
        done
    fi
fi

# ══════════════════════════════════════════════════════════════
# RESUMEN
# ══════════════════════════════════════════════════════════════
echo ""
echo "  Fase A completada"
echo "  ────────────────────────────"

n_vfdb=$(find "$ABRICATE_DIR" -name "*.vfdb.tsv" -size +0 2>/dev/null | wc -l)
n_mob=$(find "$MOBSUITE_DIR" -name "mobtyper_results.txt" 2>/dev/null | wc -l)
n_int=$(find "$INTFINDER_DIR" -name "*.integrons" 2>/dev/null | wc -l)
n_org=$(find "$AMR_ORG_DIR" -name "*_amr_org.tsv" -size +0 2>/dev/null | wc -l)
n_bakta=$(find "$BAKTA_DIR" -name "*.gff3" 2>/dev/null | wc -l)

echo "  ABRicate VFDB:       $n_vfdb ficheros"
echo "  MOB-suite:           $n_mob muestras"
echo "  IntegronFinder:      $n_int ficheros"
echo "  AMRFinderPlus --org: $n_org MAGs"
echo "  Bakta:               $n_bakta MAGs anotados"
echo ""
echo "  Para generar el Excel actualizado:"
echo "  python3 scripts/generate_excel_report.py \\"
echo "      --results-dir $RESULTS \\"
echo "      --run-name $RUN_NAME \\"
echo "      --input-dir /path/to/raw/fastqs \\"
echo "      --output $RESULTS/27_reports/${RUN_NAME}_full_report.xlsx"
