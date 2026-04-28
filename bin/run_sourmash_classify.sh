#!/bin/bash
# Clasificacion taxonomica de MAGs con Sourmash (alternativa ligera a GTDB-Tk)
# Uso: bash bin/run_sourmash_classify.sh results/260226_EPIM232 [threads]

set -euo pipefail

RESULTS="${1:?Uso: bash bin/run_sourmash_classify.sh results/<run>}"
THREADS="${2:-8}"
DB_ROOT="${DB_ROOT:-/home/asanzc/jobs/nanometa/resources}"

# Sourmash DB
SOURMASH_DB="${DB_ROOT}/sourmash/gtdb-rs226-reps-k31.zip"
LINEAGES="${DB_ROOT}/sourmash/gtdb-rs226.lineages.csv"

# Output
OUTDIR="${RESULTS}/23_taxonomy_sourmash"
mkdir -p "$OUTDIR"

echo ""
echo "  Sourmash GTDB classification"
echo "  DB: $SOURMASH_DB"
echo "  Lineages: $LINEAGES"
echo ""

# Check DB exists
if [ ! -f "$SOURMASH_DB" ]; then
    echo "ERROR: Sourmash DB no encontrada: $SOURMASH_DB"
    echo "Descargar: wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs226/gtdb-reps-rs226-k31.dna.zip"
    exit 1
fi

if [ ! -f "$LINEAGES" ]; then
    echo "ERROR: Lineages CSV no encontrado: $LINEAGES"
    exit 1
fi

for sample_dir in "$RESULTS"/21_binning_dastool/*/; do
    sample=$(basename "$sample_dir")
    bins_dir="$sample_dir/${sample}_DASTool_bins"
    [ ! -d "$bins_dir" ] && continue

    n_bins=$(ls "$bins_dir"/*.fa 2>/dev/null | wc -l)
    [ "$n_bins" -eq 0 ] && continue

    out="$OUTDIR/$sample"
    mkdir -p "$out/sketches" "$out/gather"

    echo "  $sample: $n_bins MAGs..."

    # 1. Sketch cada MAG
    for fa in "$bins_dir"/*.fa; do
        bin_name=$(basename "$fa" .fa)
        sourmash sketch dna -p k=31,scaled=1000 "$fa" \
            -o "$out/sketches/${bin_name}.sig" \
            --name "${bin_name}" 2>/dev/null || true
    done

    # 2. Gather contra GTDB
    for sig in "$out"/sketches/*.sig; do
        bin_name=$(basename "$sig" .sig)
        sourmash gather "$sig" "$SOURMASH_DB" \
            -o "$out/gather/${bin_name}.gather.csv" \
            --threshold-bp 50000 \
            -k 31 2>/dev/null || true
    done

    # 3. Concatenar gather results
    head -1 "$out"/gather/*.gather.csv 2>/dev/null | head -1 > "$out/gather_all.csv"
    tail -q -n +2 "$out"/gather/*.gather.csv >> "$out/gather_all.csv" 2>/dev/null || true

    # 4. Taxonomia
    if [ -s "$out/gather_all.csv" ]; then
        sourmash tax genome \
            -g "$out/gather_all.csv" \
            -t "$LINEAGES" \
            --output-format csv_summary \
            -o "$out/${sample}_sourmash_tax.csv" 2>/dev/null || true
    fi

    # 5. Generar TSV compatible con formato GTDB-Tk
    echo -e "user_genome\tclassification\tcontainment" > "$out/${sample}_taxonomy.tsv"
    if [ -f "$out/${sample}_sourmash_tax.csv" ]; then
        tail -n +2 "$out/${sample}_sourmash_tax.csv" | while IFS=',' read -r qname status rank fraction lineage rest; do
            echo -e "${qname}\t${lineage}\t${fraction}" >> "$out/${sample}_taxonomy.tsv"
        done
    fi

    n_classified=$(tail -n +2 "$out/${sample}_taxonomy.tsv" | wc -l)
    echo "  $sample: $n_classified MAGs clasificados"
done

echo ""
echo "  Completado. Resultados en: $OUTDIR"
echo ""
echo "  Para generar el reporte con sourmash:"
echo "  python3 scripts/generate_mag_report_v2.py \\"
echo "      --results-dir $RESULTS --run-name $(basename $RESULTS) \\"
echo "      --branding branding/ --output $RESULTS/32_reports/report.html"
