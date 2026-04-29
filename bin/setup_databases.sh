#!/bin/bash
# ╔════════════════════════════════════════════════════════════════════╗
# ║  EpiTaxMAG — Database setup helper                                 ║
# ║  Pre-populates the database root used by the Nextflow pipeline.    ║
# ╚════════════════════════════════════════════════════════════════════╝
#
# RECOMMENDED USAGE:
#
#   The pipeline downloads all required databases automatically on the first
#   run via `workflow/setup_db.nf` (`EPITAXMAG_SETUP`). Use:
#
#       nextflow run main.nf --input <fastqs> -profile <gru|censalud|pc>
#
#   This script is only needed when you want to:
#     - Pre-populate the DB root before running the pipeline
#     - Verify what is already in place (--check)
#     - Re-download a missing/corrupted file outside Nextflow
#
# USAGE:
#   bash bin/setup_databases.sh --db-root /path/to/databases       # download
#   bash bin/setup_databases.sh --db-root /path/to/databases --check  # verify only
#   bash bin/setup_databases.sh --db-root /path/to/databases --tax-only
#
# DISK FOOTPRINT (compressed before extraction):
#   TAX phase:  ~170 GB (Kraken2 78 GB + Kaiju 36 GB + Sylph ~28 GB + others)
#   MAG phase:  ~150 GB additional (GTDB-Tk 60 GB + CheckM2 1.7 GB + others)
#   The Bakta full DB (~84 GB) and AMRFinderPlus DB are downloaded by their
#   own tooling — they require the corresponding conda env to be activated.
#
# NOTES:
#   - Each DB is checked before download (skip if marker exists).
#   - wget uses -c (resumable downloads) and 5 retries.
#   - Honours http_proxy / https_proxy / NO_PROXY environment variables.
#
# ════════════════════════════════════════════════════════════════════

set -euo pipefail

# ── Colors ──────────────────────────────────────────────────────
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

# ── Arguments ──────────────────────────────────────────────────
DB_ROOT="$(dirname "$0")/../databases"
TAX_ONLY=false
CHECK_ONLY=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --db-root)  DB_ROOT="$2"; shift 2 ;;
        --tax-only) TAX_ONLY=true; shift ;;
        --check)    CHECK_ONLY=true; shift ;;
        --help|-h)
            sed -n '2,33p' "$0" | sed 's/^# \?//'
            exit 0 ;;
        *) echo -e "${RED}Unknown argument: $1${NC}" >&2; exit 1 ;;
    esac
done

DB_ROOT=$(realpath -m "$DB_ROOT")
mkdir -p "$DB_ROOT"

echo -e "${BLUE}"
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║  EpiTaxMAG — Database setup                                        ║"
echo "╠════════════════════════════════════════════════════════════════════╣"
echo "║  DB Root:   ${DB_ROOT}"
echo "║  Mode:      $(${TAX_ONLY} && echo 'TAX only (~170 GB)' || echo 'TAX + MAG')"
echo "║  Action:    $(${CHECK_ONLY} && echo 'Verify only (no downloads)' || echo 'Download missing files')"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

WGET_OPTS="-c --tries=5 --timeout=120 --show-progress"

# ── Helpers ─────────────────────────────────────────────────────

db_ok() {
    local dir="$1"
    local marker="$2"
    [[ -e "${dir}/${marker}" ]]
}

download_archive() {
    # download_archive <name> <dir> <marker> <url> <archive_basename> [strip_components]
    local name="$1" dir="$2" marker="$3" url="$4" archive="$5" strip="${6:-0}"

    if db_ok "$dir" "$marker"; then
        echo -e "  ${GREEN}✔ ${name}${NC} — already present, skipping"
        return 0
    fi

    if $CHECK_ONLY; then
        echo -e "  ${RED}✘ ${name}${NC} — missing (run without --check to download)"
        return 0
    fi

    echo -e "  ${YELLOW}↓ Downloading ${name}...${NC}"
    mkdir -p "$dir"
    cd "$dir"
    wget $WGET_OPTS -O "$archive" "$url"
    echo -e "  ${YELLOW}  Extracting ${archive}...${NC}"
    if [[ "$archive" == *.tar.gz || "$archive" == *.tgz ]]; then
        if [[ "$strip" == "0" ]]; then
            tar xzf "$archive"
        else
            tar xzf "$archive" --strip-components="$strip"
        fi
    elif [[ "$archive" == *.zip ]]; then
        unzip -q -o "$archive"
    else
        echo -e "  ${RED}Unknown archive format: $archive${NC}"
        cd - > /dev/null
        return 1
    fi
    rm -f "$archive"
    cd - > /dev/null

    if db_ok "$dir" "$marker"; then
        echo -e "  ${GREEN}✔ ${name}${NC} — installed"
    else
        echo -e "  ${RED}✘ ${name}${NC} — ERROR: marker not found after download"
        return 1
    fi
}

download_file() {
    # download_file <name> <dir> <filename> <url>
    local name="$1" dir="$2" fname="$3" url="$4"

    if db_ok "$dir" "$fname"; then
        echo -e "  ${GREEN}✔ ${name}${NC} — already present, skipping"
        return 0
    fi

    if $CHECK_ONLY; then
        echo -e "  ${RED}✘ ${name}${NC} — missing (run without --check to download)"
        return 0
    fi

    echo -e "  ${YELLOW}↓ Downloading ${name}...${NC}"
    mkdir -p "$dir"
    wget $WGET_OPTS -O "${dir}/${fname}" "$url"

    if db_ok "$dir" "$fname"; then
        echo -e "  ${GREEN}✔ ${name}${NC} — installed"
    else
        echo -e "  ${RED}✘ ${name}${NC} — ERROR: download failed"
        return 1
    fi
}

# ══════════════════════════════════════════════════════════════════
# TAX-phase databases
# ══════════════════════════════════════════════════════════════════

echo -e "\n${BLUE}── TAX phase ────────────────────────────────────────${NC}\n"

# 1. Kraken2 PlusPF (~78 GB compressed)
download_archive "Kraken2 PlusPF" \
    "${DB_ROOT}/k2_pluspf_20251015" \
    "hash.k2d" \
    "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20251015.tar.gz" \
    "k2_pluspf_20251015.tar.gz" || true

# 2. Kaiju RefSeq ref 2024-08-14 (~36 GB compressed)
download_archive "Kaiju RefSeq ref" \
    "${DB_ROOT}/kaiju" \
    "kaiju_db_refseq_ref.fmi" \
    "https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_refseq_ref_2024-08-14.tgz" \
    "kaiju_db_refseq_ref_2024-08-14.tgz" || true

# 3. Sylph databases (CMU primary, Google Cloud mirror)
download_file "Sylph GTDB r226" \
    "${DB_ROOT}/sylph" \
    "gtdb-r226-c200-dbv1.syldb" \
    "http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r226-c200-dbv1.syldb" || true

download_file "Sylph Fungi 2025-10-11" \
    "${DB_ROOT}/sylph" \
    "fungi-refseq-2025-10-11-c200-dbv1.syldb" \
    "http://faust.compbio.cs.cmu.edu/sylph-stuff/fungi-refseq-2025-10-11-c200-dbv1.syldb" || true

download_file "Sylph Viral IMG/VR" \
    "${DB_ROOT}/sylph" \
    "imgvr_c200_v0.3.0.syldb" \
    "http://faust.compbio.cs.cmu.edu/sylph-stuff/imgvr_c200_v0.3.0.syldb" || true

# 4. ResFinder DB (~3 MB) — git clone, then concatenate
RESFINDER_DIR="${DB_ROOT}/kma_resfinder/resfinder_db"
if db_ok "$RESFINDER_DIR" "all.fsa"; then
    echo -e "  ${GREEN}✔ ResFinder${NC} — already present, skipping"
elif $CHECK_ONLY; then
    echo -e "  ${RED}✘ ResFinder${NC} — missing (run without --check to download)"
else
    echo -e "  ${YELLOW}↓ Cloning ResFinder DB...${NC}"
    mkdir -p "$RESFINDER_DIR"
    rm -rf "${RESFINDER_DIR}_tmp"
    git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git "${RESFINDER_DIR}_tmp"
    cat "${RESFINDER_DIR}_tmp"/*.fsa > "${RESFINDER_DIR}/all.fsa"
    cp "${RESFINDER_DIR}_tmp/phenotypes.txt" "$RESFINDER_DIR/" 2>/dev/null || true
    cp "${RESFINDER_DIR}_tmp/config"         "$RESFINDER_DIR/" 2>/dev/null || true
    rm -rf "${RESFINDER_DIR}_tmp"
    echo -e "  ${GREEN}✔ ResFinder${NC} — installed"
fi

# 5. FastQ Screen panel — only generate a placeholder config; the panel itself
# is hosted on Zenodo and is fetched by the Nextflow setup workflow.
FQS_DIR="${DB_ROOT}/fastqscreen"
if db_ok "$FQS_DIR" "fastq_screen.conf"; then
    echo -e "  ${GREEN}✔ FastQ Screen${NC} — already configured"
else
    echo -e "  ${YELLOW}⚠ FastQ Screen${NC} — not configured here."
    echo -e "    The panel is downloaded by 'nextflow run main.nf' from:"
    echo -e "    https://zenodo.org/records/19860716/files/epitaxmag-fqs-panel-v1.tar.gz"
    echo -e "    To configure manually, place .fasta files under: ${FQS_DIR}/genomes/"
fi

# ══════════════════════════════════════════════════════════════════
# MAG-phase databases
# ══════════════════════════════════════════════════════════════════

if ! $TAX_ONLY; then

echo -e "\n${BLUE}── MAG phase ────────────────────────────────────────${NC}\n"

# 6. CheckM2 (~1.7 GB)
download_archive "CheckM2" \
    "${DB_ROOT}/checkm2" \
    "CheckM2_database/uniref100.KO.1.dmnd" \
    "https://zenodo.org/records/5571251/files/checkm2_database.tar.gz" \
    "checkm2_database.tar.gz" || true

# 7. GTDB-Tk r232 (~60 GB compressed)
download_archive "GTDB-Tk r232" \
    "${DB_ROOT}/gtdbtk_r232" \
    "release232" \
    "https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz" \
    "gtdbtk_r232_data.tar.gz" || true

# 8. geNomad v1.9 (~840 MB)
download_archive "geNomad v1.9" \
    "${DB_ROOT}/genomad" \
    "genomad_db" \
    "https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz" \
    "genomad_db_v1.9.tar.gz" || true

# 9. AMRFinderPlus and Bakta require their own conda envs:
echo ""
echo -e "  ${YELLOW}AMRFinderPlus and Bakta DBs require their own tooling:${NC}"
echo -e "    The Nextflow setup workflow downloads them automatically inside"
echo -e "    the matching conda env. To pre-fetch manually, activate the env"
echo -e "    and run:"
echo -e "      conda activate envs/nf-mag-amrfinderplus.yml"
echo -e "      amrfinder_update --database ${DB_ROOT}/amrfinderplus/latest"
echo -e ""
echo -e "      conda activate envs/nf-mag-bakta.yml"
echo -e "      bakta_db download --output ${DB_ROOT}/bakta/db --type full"

fi  # MAG phase

# ══════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════

echo -e "\n${BLUE}── Summary ──────────────────────────────────────────${NC}\n"
echo -e "  DB Root: ${DB_ROOT}"
echo ""

TOTAL_OK=0
TOTAL_MISSING=0

check_and_count() {
    local name="$1" dir="$2" marker="$3"
    if db_ok "$dir" "$marker"; then
        echo -e "  ${GREEN}✔${NC} $name"
        TOTAL_OK=$((TOTAL_OK + 1))
    else
        echo -e "  ${RED}✘${NC} $name"
        TOTAL_MISSING=$((TOTAL_MISSING + 1))
    fi
}

echo "  TAX phase:"
check_and_count "Kraken2 PlusPF"     "${DB_ROOT}/k2_pluspf_20251015"           "hash.k2d"
check_and_count "Kaiju RefSeq ref"   "${DB_ROOT}/kaiju"                        "kaiju_db_refseq_ref.fmi"
check_and_count "Sylph GTDB"         "${DB_ROOT}/sylph"                        "gtdb-r226-c200-dbv1.syldb"
check_and_count "Sylph Fungi"        "${DB_ROOT}/sylph"                        "fungi-refseq-2025-10-11-c200-dbv1.syldb"
check_and_count "Sylph Viral"        "${DB_ROOT}/sylph"                        "imgvr_c200_v0.3.0.syldb"
check_and_count "ResFinder"          "${DB_ROOT}/kma_resfinder/resfinder_db"   "all.fsa"
check_and_count "FastQ Screen"       "${DB_ROOT}/fastqscreen"                  "fastq_screen.conf"

if ! $TAX_ONLY; then
echo "  MAG phase:"
check_and_count "GTDB-Tk r232"       "${DB_ROOT}/gtdbtk_r232"                  "release232"
check_and_count "CheckM2"            "${DB_ROOT}/checkm2"                      "CheckM2_database/uniref100.KO.1.dmnd"
check_and_count "geNomad"            "${DB_ROOT}/genomad"                      "genomad_db"
check_and_count "AMRFinderPlus"      "${DB_ROOT}/amrfinderplus"                "latest"
check_and_count "Bakta"              "${DB_ROOT}/bakta"                        "db"
fi

echo ""
echo -e "  ${GREEN}${TOTAL_OK} OK${NC} | ${RED}${TOTAL_MISSING} missing${NC}"

if [[ $TOTAL_MISSING -gt 0 ]]; then
    echo ""
    echo -e "  ${YELLOW}Run the pipeline (which auto-downloads any missing DB):${NC}"
    echo -e "    nextflow run main.nf --input data/... --db_root ${DB_ROOT} -profile gru"
else
    echo ""
    echo -e "  ${GREEN}All databases ready.${NC}"
    echo -e "    nextflow run main.nf --input data/... --db_root ${DB_ROOT} -profile gru"
fi
echo ""
