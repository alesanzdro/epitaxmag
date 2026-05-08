#!/bin/bash
# Build KMA database from type-strain genomes for phenotypic verification
# Usage: bash bin/build_phenotypic_kma_db.sh

set -euo pipefail

OUTDIR="assets/phenotypic_targets"
mkdir -p "$OUTDIR"

# Proxy honoured from the environment if set; this script does not
# hard-code institutional proxies. Configure http_proxy/https_proxy in
# your shell or in conf/local.config (env.http_proxy, env.https_proxy).

echo "  Building phenotypic type-strain KMA database"
echo ""

# Type-strain accessions (RefSeq reference genomes)
declare -A TARGETS=(
    ["Escherichia_coli"]="GCF_000005845.2"
    ["Klebsiella_pneumoniae"]="GCF_000240185.1"
    ["Klebsiella_oxytoca"]="GCF_003261575.2"
    ["Enterobacter_cloacae"]="GCF_000025565.1"
    ["Citrobacter_freundii"]="GCF_003697165.2"
    ["Pseudomonas_aeruginosa"]="GCF_000006765.1"
    ["Klebsiella_aerogenes"]="GCF_007632255.1"
    ["Enterococcus_faecium"]="GCF_000174395.2"
    ["Enterococcus_faecalis"]="GCF_000007785.1"
    ["Salmonella_enterica"]="GCF_000006945.2"
)

rm -f "$OUTDIR/phenotypic_targets.fa"

for species in "${!TARGETS[@]}"; do
    acc="${TARGETS[$species]}"
    echo "  Downloading $species ($acc)..."

    # Download using NCBI datasets API
    wget -q "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${acc}/download?include_annotation_type=GENOME_FASTA" \
        -O "$OUTDIR/${species}.zip" 2>/dev/null || \
    wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/$(echo $acc | sed 's/GCF_/GCF\//' | sed 's/\(.\{3\}\)/\1\//g' | sed 's/\/$//')/${acc}*/${acc}*_genomic.fna.gz" \
        -O "$OUTDIR/${species}.fna.gz" 2>/dev/null || true

    # Try unzip (datasets format)
    if [ -f "$OUTDIR/${species}.zip" ] && [ -s "$OUTDIR/${species}.zip" ]; then
        unzip -qjo "$OUTDIR/${species}.zip" "*.fna" -d "$OUTDIR/tmp_${species}/" 2>/dev/null || true
        rm -f "$OUTDIR/${species}.zip"
    elif [ -f "$OUTDIR/${species}.fna.gz" ] && [ -s "$OUTDIR/${species}.fna.gz" ]; then
        mkdir -p "$OUTDIR/tmp_${species}"
        gunzip -c "$OUTDIR/${species}.fna.gz" > "$OUTDIR/tmp_${species}/${species}.fna"
        rm -f "$OUTDIR/${species}.fna.gz"
    fi

    # Rewrite headers with species name
    if ls "$OUTDIR/tmp_${species}/"*.fna 1>/dev/null 2>&1; then
        awk -v sp="${species}" '/^>/ { sub(/^>/, ">" sp " "); print; next } { print }' \
            "$OUTDIR/tmp_${species}/"*.fna >> "$OUTDIR/phenotypic_targets.fa"
        echo "  OK: $species added"
    else
        echo "  WARN: Could not download $species"
    fi
    rm -rf "$OUTDIR/tmp_${species}"
done

# Index with KMA
if [ -f "$OUTDIR/phenotypic_targets.fa" ] && [ -s "$OUTDIR/phenotypic_targets.fa" ]; then
    echo ""
    echo "  Indexing with KMA..."
    kma index -i "$OUTDIR/phenotypic_targets.fa" -o "$OUTDIR/phenotypic_targets_db" 2>/dev/null || \
    echo "  WARN: kma not in PATH. Index manually: kma index -i $OUTDIR/phenotypic_targets.fa -o $OUTDIR/phenotypic_targets_db"

    echo ""
    n_species=$(grep -c "^>" "$OUTDIR/phenotypic_targets.fa")
    echo "  Database built: $n_species contigs from ${#TARGETS[@]} species"
    ls -lh "$OUTDIR/phenotypic_targets"*
else
    echo "  ERROR: No genomes downloaded. Check network/proxy."
fi
