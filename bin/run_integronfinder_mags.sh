#!/bin/bash
# IntegronFinder over MAGs (bins) — much faster than over full assemblies
# Usage: bash bin/run_integronfinder_mags.sh results/260226_EPIM232

set -euo pipefail

RESULTS="${1:?Usage: bash bin/run_integronfinder_mags.sh results/<run>}"
THREADS="${2:-8}"
MAX_PARALLEL="${3:-6}"

CACHE="${NXF_SINGULARITY_CACHEDIR:-/GRU0/DATABASES/singularity_cache}"
BIND="--bind /GRU0 --bind /home/asanzc/jobs"
SIF="$CACHE/depot.galaxyproject.org-singularity-integron_finder%3A2.0rc6--py_0.img"
OUTDIR="$RESULTS/32_integronfinder_mags"
TMPBASE="$RESULTS/../work/tmp"

mkdir -p "$OUTDIR" "$TMPBASE"

echo "  IntegronFinder over MAGs: $MAX_PARALLEL parallel x $THREADS threads"

RUNNING=0
TOTAL=0
for sample_dir in "$RESULTS"/21_binning_dastool/*/; do
    sample=$(basename "$sample_dir")
    bins_dir="$sample_dir/${sample}_DASTool_bins"
    [ ! -d "$bins_dir" ] && continue

    for bin_fa in "$bins_dir"/*.fa; do
        [ ! -f "$bin_fa" ] && continue
        bin_name=$(basename "$bin_fa" .fa)
        out_dir="$OUTDIR/${sample}/${bin_name}"

        # Skip if already complete
        if [ -d "$out_dir" ] && find "$out_dir" -name "*.integrons" 2>/dev/null | grep -q .; then
            continue
        fi

        rm -rf "$out_dir"
        mkdir -p "$out_dir"
        TOTAL=$((TOTAL + 1))

        (
            export TMPDIR="$TMPBASE/intfinder_${bin_name}"
            mkdir -p "$TMPDIR"
            singularity exec --no-home $BIND "$SIF" \
                integron_finder \
                    --local-max --func-annot \
                    --cpu "$THREADS" \
                    --outdir "$out_dir" \
                    "$bin_fa" >"$out_dir/run.log" 2>&1
            rm -rf "$TMPDIR"
        ) &

        RUNNING=$((RUNNING + 1))
        if [ "$RUNNING" -ge "$MAX_PARALLEL" ]; then
            wait -n
            RUNNING=$((RUNNING - 1))
        fi
    done
done

wait

echo ""
echo "  Done: $TOTAL MAGs processed"
n_found=$(find "$OUTDIR" -name "*.integrons" -size +0 2>/dev/null | xargs grep -l "^[^#]" 2>/dev/null | wc -l)
echo "  MAGs with integrons: $n_found"
