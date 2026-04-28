#!/bin/bash
# EpiTaxMAG Fase A — Modulos pendientes: AMRFinderPlus --organism + Bakta
# Uso: bash bin/run_phase_a_remaining.sh results/260226_EPIM232 16

set -euo pipefail

RESULTS="${1:?Uso: bash bin/run_phase_a_remaining.sh results/<run>}"
THREADS="${2:-16}"

CACHE="${NXF_SINGULARITY_CACHEDIR:-/GRU0/DATABASES/singularity_cache}"
BIND="--bind /GRU0 --bind /home/asanzc/jobs"
DB_ROOT="${DB_ROOT:-/home/asanzc/jobs/nanometa/resources}"

SIF_AMRFINDER="$CACHE/depot.galaxyproject.org-singularity-ncbi-amrfinderplus%3A4.2.7--hf69ffd2_0.img"
SIF_BAKTA="$CACHE/depot.galaxyproject.org-singularity-bakta%3A1.12.0--pyhdfd78af_0.img"

FAKEHOME="$(realpath "$RESULTS/../work/tmp/fakehome")"
mkdir -p "$FAKEHOME"
srun() { singularity exec --home "$FAKEHOME" $BIND "$@"; }

echo ""
echo "  Fase A restante: AMRFinderPlus --organism + Bakta"
echo "  Threads: $THREADS"
echo ""

# ══════════════════════════════════════════════════════════════
# MODULO 1: AMRFinderPlus --organism (mutaciones puntuales)
# ══════════════════════════════════════════════════════════════
AMR_ORG_DIR="$RESULTS/33_amr_organism"
echo "[1/2] AMRFinderPlus --organism..."

mkdir -p "$AMR_ORG_DIR"

genus_to_organism() {
    case "$1" in
        Acinetobacter)   echo "Acinetobacter_baumannii" ;;
        Escherichia)     echo "Escherichia" ;;
        Klebsiella)      echo "Klebsiella_pneumoniae" ;;
        Pseudomonas)     echo "Pseudomonas_aeruginosa" ;;
        Salmonella)      echo "Salmonella" ;;
        Staphylococcus)  echo "Staphylococcus_aureus" ;;
        Enterococcus)    echo "Enterococcus_faecalis" ;;
        Neisseria)       echo "Neisseria_meningitidis" ;;
        Campylobacter)   echo "Campylobacter" ;;
        Clostridioides)  echo "Clostridioides_difficile" ;;
        Burkholderia)    echo "Burkholderia_cepacia" ;;
        Vibrio)          echo "Vibrio_cholerae" ;;
        *)               echo "" ;;
    esac
}

for sample_dir in "$RESULTS"/21_binning_dastool/*/; do
    sample=$(basename "$sample_dir")
    bins_dir="$sample_dir/${sample}_DASTool_bins"
    tax_file="$RESULTS/23_taxonomy_gtdbtk/$sample/${sample}_taxonomy.tsv"

    [ ! -d "$bins_dir" ] && continue
    [ ! -f "$tax_file" ] && continue

    mkdir -p "$AMR_ORG_DIR/$sample"

    while IFS=$'\t' read -r bin_name taxonomy rest; do
        [ "$bin_name" = "user_genome" ] && continue
        bin_fa="$bins_dir/${bin_name}.fa"
        [ ! -f "$bin_fa" ] && continue

        out="$AMR_ORG_DIR/$sample/${bin_name}_amr_org.tsv"
        [ -f "$out" ] && [ -s "$out" ] && continue

        genus=$(echo "$taxonomy" | grep -oP 'g__\K[^;]+' || echo "")
        organism=$(genus_to_organism "$genus")

        if [ -n "$organism" ]; then
            echo "  $bin_name -> $organism"
            srun "$SIF_AMRFINDER" amrfinder \
                --nucleotide "$bin_fa" \
                --organism "$organism" \
                --database "$DB_ROOT/amrfinderplus/latest" \
                --threads "$THREADS" --plus \
                --name "$bin_name" \
                --output "$out" 2>/dev/null || true
        else
            srun "$SIF_AMRFINDER" amrfinder \
                --nucleotide "$bin_fa" \
                --database "$DB_ROOT/amrfinderplus/latest" \
                --threads "$THREADS" --plus \
                --name "$bin_name" \
                --output "$out" 2>/dev/null || true
        fi
    done < "$tax_file"
    echo "  $sample: completado"
done
echo "  AMRFinderPlus --organism: OK"

# ══════════════════════════════════════════════════════════════
# MODULO 2: Bakta — Anotacion funcional de MAGs
# ══════════════════════════════════════════════════════════════
BAKTA_DIR="$RESULTS/34_bakta"
BAKTA_DB="$DB_ROOT/bakta/db"
echo ""
echo "[2/2] Bakta (anotacion funcional)..."

if [ ! -d "$BAKTA_DB" ]; then
    echo "  ERROR: Bakta DB no encontrada en $BAKTA_DB"
    exit 1
fi

mkdir -p "$BAKTA_DIR"
export TMPDIR="$RESULTS/../work/tmp"
mkdir -p "$TMPDIR"

for sample_dir in "$RESULTS"/21_binning_dastool/*/; do
    sample=$(basename "$sample_dir")
    bins_dir="$sample_dir/${sample}_DASTool_bins"
    [ ! -d "$bins_dir" ] && continue
    n_bins=$(ls "$bins_dir"/*.fa 2>/dev/null | wc -l)
    [ "$n_bins" -eq 0 ] && continue

    mkdir -p "$BAKTA_DIR/$sample"
    DONE_BAKTA=0

    for bin_fa in "$bins_dir"/*.fa; do
        bin_name=$(basename "$bin_fa" .fa)
        out_dir="$BAKTA_DIR/$sample/$bin_name"

        # Skip si ya completado
        if [ -d "$out_dir" ] && [ -f "$out_dir/${bin_name}.gff3" ]; then
            DONE_BAKTA=$((DONE_BAKTA+1))
            continue
        fi

        echo "  $bin_name..."
        mkdir -p "$out_dir"
        srun "$SIF_BAKTA" bash -c "
            export TMPDIR=/home/asanzc/jobs/epitaxmag/work/tmp
            bakta \
                --db $BAKTA_DB \
                --output $out_dir \
                --prefix $bin_name \
                --threads $THREADS \
                --locus-tag ${sample} \
                --skip-plot \
                --skip-cds \
                --force \
                $bin_fa" 2>"$out_dir/bakta.log" || true
        DONE_BAKTA=$((DONE_BAKTA+1))
    done
    echo "  $sample: $DONE_BAKTA/$n_bins MAGs anotados"
done

# ══════════════════════════════════════════════════════════════
# RESUMEN
# ══════════════════════════════════════════════════════════════
echo ""
echo "  Fase A restante completada"
echo "  AMRFinderPlus --organism: $(find "$AMR_ORG_DIR" -name "*_amr_org.tsv" -size +0 2>/dev/null | wc -l) MAGs"
echo "  Bakta: $(find "$BAKTA_DIR" -name "*.gff3" 2>/dev/null | wc -l) MAGs anotados"
echo ""
echo "  Regenerar Excel:"
echo "  python3 scripts/generate_excel_report.py --results-dir $RESULTS --run-name $(basename $RESULTS) --input-dir /path/to/raw/fastqs --output $RESULTS/27_reports/$(basename $RESULTS)_full_report_final.xlsx"
