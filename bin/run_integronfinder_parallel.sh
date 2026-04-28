#!/bin/bash
# IntegronFinder en paralelo — 4 muestras simultaneas
# Uso: bash bin/run_integronfinder_parallel.sh results/260226_EPIM232

set -euo pipefail

RESULTS="${1:?Uso: bash bin/run_integronfinder_parallel.sh results/<run>}"
THREADS_PER_JOB="${2:-12}"
MAX_PARALLEL="${3:-4}"

CACHE="${NXF_SINGULARITY_CACHEDIR:-/GRU0/DATABASES/singularity_cache}"
BIND="--bind /GRU0 --bind /home/asanzc/jobs"
SIF="$CACHE/depot.galaxyproject.org-singularity-integron_finder%3A2.0rc6--py_0.img"
OUTDIR="$RESULTS/32_integronfinder"
TMPBASE="$RESULTS/../work/tmp"

mkdir -p "$OUTDIR" "$TMPBASE"

echo "  IntegronFinder paralelo: $MAX_PARALLEL jobs x $THREADS_PER_JOB threads"
echo ""

RUNNING=0
for asm in "$RESULTS"/15_polish_medaka/*.polished.fasta; do
    [ ! -f "$asm" ] && continue
    sample=$(basename "$asm" .polished.fasta)
    out_dir="$OUTDIR/$sample"

    # Skip si ya tiene resultados completos
    if [ -d "$out_dir" ] && find "$out_dir" -name "*.integrons" -size +0 2>/dev/null | grep -q .; then
        echo "  $sample: ya completado, saltando"
        continue
    fi

    # Limpiar resultados parciales
    rm -rf "$out_dir"
    mkdir -p "$out_dir"

    echo "  $sample: lanzando IntegronFinder..."
    (
        export TMPDIR="$TMPBASE/intfinder_${sample}"
        mkdir -p "$TMPDIR"
        singularity exec --no-home $BIND "$SIF" \
            integron_finder \
                --local-max --func-annot \
                --cpu "$THREADS_PER_JOB" \
                --outdir "$out_dir" \
                "$asm" 2>"$out_dir/integronfinder.log"
        rm -rf "$TMPDIR"
        echo "  $sample: COMPLETADO"
    ) &

    RUNNING=$((RUNNING + 1))
    if [ "$RUNNING" -ge "$MAX_PARALLEL" ]; then
        wait -n
        RUNNING=$((RUNNING - 1))
    fi
done

wait
echo ""
echo "  IntegronFinder completado para todas las muestras"
echo "  Resultados en: $OUTDIR"
