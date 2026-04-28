#!/bin/bash
# ╔════════════════════════════════════════════════════════════════════╗
# ║  EpiTaxMAG — Setup de Bases de Datos                             ║
# ║  Descarga y configura todas las DBs necesarias para el pipeline  ║
# ╚════════════════════════════════════════════════════════════════════╝
#
# USO:
#   bash bin/setup_databases.sh                          # Default: ./databases
#   bash bin/setup_databases.sh --db-root /disco/dbs     # Ruta personalizada
#   bash bin/setup_databases.sh --tax-only               # Solo fase TAX (~170 GB)
#   bash bin/setup_databases.sh --check                  # Verificar sin descargar
#
# ESPACIO NECESARIO:
#   Fase TAX:  ~170 GB (Kraken2 100GB + Kaiju 34GB + Sylph 17GB + otros)
#   Fase MAG:  ~600 GB (GTDB-Tk 208GB + Bakta 84GB + eggNOG 48GB + otros)
#   Total:     ~780 GB
#
# NOTAS:
#   - Cada DB se comprueba antes de descargar (skip si existe)
#   - wget con -c (reanuda descargas interrumpidas)
#   - Compatible con proxies (configura http_proxy/https_proxy)
#   - El pipeline detecta automáticamente las DBs en --db_root
#
# ════════════════════════════════════════════════════════════════════

set -euo pipefail

# ── Colores ─────────────────────────────────────────────────────
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

# ── Argumentos ──────────────────────────────────────────────────
DB_ROOT="$(dirname "$0")/../databases"
TAX_ONLY=false
CHECK_ONLY=false
THREADS=8

while [[ $# -gt 0 ]]; do
    case "$1" in
        --db-root)  DB_ROOT="$2"; shift 2 ;;
        --tax-only) TAX_ONLY=true; shift ;;
        --check)    CHECK_ONLY=true; shift ;;
        --threads)  THREADS="$2"; shift 2 ;;
        --help|-h)
            head -25 "$0" | grep "^#" | sed 's/^# \?//'
            exit 0 ;;
        *) echo -e "${RED}Argumento desconocido: $1${NC}"; exit 1 ;;
    esac
done

DB_ROOT=$(realpath -m "$DB_ROOT")
mkdir -p "$DB_ROOT"

echo -e "${BLUE}"
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║  EpiTaxMAG — Setup de Bases de Datos                             ║"
echo "╠════════════════════════════════════════════════════════════════════╣"
echo "║  DB Root:   ${DB_ROOT}"
echo "║  Modo:      $(${TAX_ONLY} && echo 'Solo TAX (~170 GB)' || echo 'Completo TAX+MAG (~780 GB)')"
echo "║  Check:     $(${CHECK_ONLY} && echo 'Solo verificación' || echo 'Descarga si falta')"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

WGET_OPTS="-c --show-progress --timeout=60 --tries=3"

# ── Funciones ───────────────────────────────────────────────────

db_ok() {
    local dir="$1"
    local marker="$2"
    [[ -e "${dir}/${marker}" ]]
}

status() {
    local name="$1"
    local dir="$2"
    local marker="$3"
    if db_ok "$dir" "$marker"; then
        echo -e "  ${GREEN}✔ ${name}${NC} — listo"
        return 0
    else
        echo -e "  ${RED}✘ ${name}${NC} — no encontrado"
        return 1
    fi
}

download_if_missing() {
    local name="$1"
    local dir="$2"
    local marker="$3"
    local url="$4"
    local extract_cmd="${5:-}"

    if db_ok "$dir" "$marker"; then
        echo -e "  ${GREEN}✔ ${name}${NC} — ya descargado, saltando"
        return 0
    fi

    if $CHECK_ONLY; then
        echo -e "  ${RED}✘ ${name}${NC} — falta (ejecutar sin --check para descargar)"
        return 1
    fi

    echo -e "  ${YELLOW}↓ Descargando ${name}...${NC}"
    mkdir -p "$dir"

    local filename=$(basename "$url")
    wget $WGET_OPTS -P "$dir" "$url"

    if [[ -n "$extract_cmd" ]]; then
        echo -e "  ${YELLOW}  Extrayendo...${NC}"
        eval "$extract_cmd"
    fi

    if db_ok "$dir" "$marker"; then
        echo -e "  ${GREEN}✔ ${name}${NC} — instalado correctamente"
    else
        echo -e "  ${RED}✘ ${name}${NC} — ERROR: marcador no encontrado tras descarga"
        return 1
    fi
}

# ══════════════════════════════════════════════════════════════════
# FASE TAX — Bases de datos para perfilado taxonómico
# ══════════════════════════════════════════════════════════════════

echo -e "\n${BLUE}── FASE TAX ──────────────────────────────────────────${NC}\n"

# 1. Kraken2 PlusPF (~100 GB)
KRAKEN2_DIR="${DB_ROOT}/k2_pluspf_20251015"
download_if_missing \
    "Kraken2 PlusPF" \
    "$KRAKEN2_DIR" \
    "hash.k2d" \
    "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20251015.tar.gz" \
    "tar xzf ${KRAKEN2_DIR}/k2_pluspf_20251015.tar.gz -C ${KRAKEN2_DIR}/ && rm -f ${KRAKEN2_DIR}/k2_pluspf_20251015.tar.gz" \
    || true

# 2. Kaiju RefSeq (~34 GB)
KAIJU_DIR="${DB_ROOT}/kaiju"
download_if_missing \
    "Kaiju RefSeq" \
    "$KAIJU_DIR" \
    "kaiju_db_refseq_ref.fmi" \
    "https://bioinformatics-centre.github.io/kaiju/downloads/kaiju_db_refseq_ref.tgz" \
    "tar xzf ${KAIJU_DIR}/kaiju_db_refseq_ref.tgz -C ${KAIJU_DIR}/ && rm -f ${KAIJU_DIR}/kaiju_db_refseq_ref.tgz" \
    || true

# 3. Sylph GTDB r226 (~14 GB)
SYLPH_DIR="${DB_ROOT}/sylph"
download_if_missing \
    "Sylph GTDB r226" \
    "$SYLPH_DIR" \
    "gtdb-r226-c200-dbv1.syldb" \
    "https://storage.googleapis.com/sylph-stuff/v0.7-databases/gtdb-r226-c200-dbv1.syldb" \
    || true

# 4. Sylph Fungi RefSeq (~700 MB)
download_if_missing \
    "Sylph Fungi" \
    "$SYLPH_DIR" \
    "fungi-refseq-2024-07-25-c200-v0.3.syldb" \
    "https://storage.googleapis.com/sylph-stuff/v0.7-databases/fungi-refseq-2024-07-25-c200-v0.3.syldb" \
    || true

# 5. Sylph IMG/VR Viral (~2 GB)
download_if_missing \
    "Sylph Viral IMG/VR" \
    "$SYLPH_DIR" \
    "imgvr_c200_v0.3.0.syldb" \
    "https://storage.googleapis.com/sylph-stuff/v0.7-databases/imgvr_c200_v0.3.0.syldb" \
    || true

# 6. ResFinder/KMA (~3 MB)
RESFINDER_DIR="${DB_ROOT}/kma_resfinder/resfinder_db"
if ! db_ok "$RESFINDER_DIR" "all.fsa"; then
    if ! $CHECK_ONLY; then
        echo -e "  ${YELLOW}↓ Descargando ResFinder DB...${NC}"
        mkdir -p "$RESFINDER_DIR"
        git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git "${RESFINDER_DIR}_tmp" 2>/dev/null || true
        if [[ -d "${RESFINDER_DIR}_tmp" ]]; then
            # Concatenar todas las secuencias en all.fsa
            cat "${RESFINDER_DIR}_tmp"/*.fsa > "${RESFINDER_DIR}/all.fsa" 2>/dev/null || true
            cp "${RESFINDER_DIR}_tmp/phenotypes.txt" "$RESFINDER_DIR/" 2>/dev/null || true
            cp "${RESFINDER_DIR}_tmp/config" "$RESFINDER_DIR/" 2>/dev/null || true
            rm -rf "${RESFINDER_DIR}_tmp"
            echo -e "  ${GREEN}✔ ResFinder${NC} — instalado"
        fi
    else
        echo -e "  ${RED}✘ ResFinder${NC} — falta"
    fi
else
    echo -e "  ${GREEN}✔ ResFinder${NC} — ya descargado, saltando"
fi

# 7. FastQ Screen genomes (requiere configuración manual)
FQS_DIR="${DB_ROOT}/FastqScreen"
if ! db_ok "$FQS_DIR" "fastq_screen.conf"; then
    echo -e "  ${YELLOW}⚠ FastQ Screen${NC} — requiere configuración manual"
    echo -e "    Crear: ${FQS_DIR}/fastq_screen.conf"
    echo -e "    Con genomas de referencia en: ${FQS_DIR}/genomes/"
    echo -e "    Ver: docs/epitax_pipeline.md para detalles"
    if ! $CHECK_ONLY; then
        mkdir -p "${FQS_DIR}/genomes"
        # Crear config mínimo con H. sapiens si no existe
        if [[ ! -f "${FQS_DIR}/fastq_screen.conf" ]]; then
            cat > "${FQS_DIR}/fastq_screen.conf" << 'CONF'
# FastQ Screen configuration — EpiTaxMAG
# Editar las rutas DATABASE según tus genomas de referencia
ALIGNER   minimap2
THREADS   8

# Ejemplo (descomentar y ajustar rutas):
# DATABASE  H_sapiens  /path/to/genomes/H_sapiens
# DATABASE  Ecoli      /path/to/genomes/Ecoli
CONF
            echo -e "    Config plantilla creado en ${FQS_DIR}/fastq_screen.conf"
        fi
    fi
else
    echo -e "  ${GREEN}✔ FastQ Screen${NC} — configurado"
fi

# ══════════════════════════════════════════════════════════════════
# FASE MAG — Bases de datos para ensamblaje y MAGs
# ══════════════════════════════════════════════════════════════════

if ! $TAX_ONLY; then

echo -e "\n${BLUE}── FASE MAG ──────────────────────────────────────────${NC}\n"

# 8. GTDB-Tk r220 (~208 GB)
GTDBTK_DIR="${DB_ROOT}/gtdbtk"
download_if_missing \
    "GTDB-Tk r220" \
    "$GTDBTK_DIR" \
    "release220/metadata/metadata.txt" \
    "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz" \
    "tar xzf ${GTDBTK_DIR}/gtdbtk_r220_data.tar.gz -C ${GTDBTK_DIR}/ && rm -f ${GTDBTK_DIR}/gtdbtk_r220_data.tar.gz" \
    || true

# 9. CheckM2 (~3 GB)
CHECKM2_DIR="${DB_ROOT}/checkm2"
download_if_missing \
    "CheckM2" \
    "$CHECKM2_DIR" \
    "CheckM2_database/uniref100.KO.1.dmnd" \
    "https://zenodo.org/records/5571251/files/checkm2_database.tar.gz" \
    "tar xzf ${CHECKM2_DIR}/checkm2_database.tar.gz -C ${CHECKM2_DIR}/ && rm -f ${CHECKM2_DIR}/checkm2_database.tar.gz" \
    || true

# 10. geNomad (~1.4 GB)
GENOMAD_DIR="${DB_ROOT}/genomad"
download_if_missing \
    "geNomad" \
    "$GENOMAD_DIR" \
    "genomad_db/genomad_db" \
    "https://zenodo.org/records/8339387/files/genomad_db_v1.7.tar.gz" \
    "tar xzf ${GENOMAD_DIR}/genomad_db_v1.7.tar.gz -C ${GENOMAD_DIR}/ && rm -f ${GENOMAD_DIR}/genomad_db_v1.7.tar.gz" \
    || true

# 11. CheckV (~6.4 GB)
CHECKV_DIR="${DB_ROOT}/checkv"
download_if_missing \
    "CheckV" \
    "$CHECKV_DIR" \
    "checkv-db-v1.5/README.txt" \
    "https://zenodo.org/records/7899058/files/checkv-db-v1.5.tar.gz" \
    "tar xzf ${CHECKV_DIR}/checkv-db-v1.5.tar.gz -C ${CHECKV_DIR}/ && rm -f ${CHECKV_DIR}/checkv-db-v1.5.tar.gz" \
    || true

# Herramientas que necesitan sus propios comandos de descarga:
echo -e "\n  ${YELLOW}Las siguientes DBs requieren sus herramientas para descargar:${NC}"
echo -e "  Si usas Singularity, activa los containers correspondientes."
echo ""

for db_info in \
    "Bakta:bakta:bakta db download --output ${DB_ROOT}/bakta/db --type full" \
    "BUSCO:busco:busco --download prokaryota --download_path ${DB_ROOT}/busco" \
    "GUNC:gunc:gunc download_db ${DB_ROOT}/gunc" \
    "eggNOG:eggnog-mapper:download_eggnog_data.py -y --data_dir ${DB_ROOT}/eggnog" \
    "AMRFinderPlus:amrfinderplus:amrfinder_update --database ${DB_ROOT}/amrfinderplus/latest"
do
    IFS=':' read -r name tool cmd <<< "$db_info"
    dir="${DB_ROOT}/$(echo "$name" | tr '[:upper:]' '[:lower:]')"
    if [[ -d "$dir" ]] && [[ $(find "$dir" -type f 2>/dev/null | head -1) ]]; then
        echo -e "  ${GREEN}✔ ${name}${NC} — ya existe"
    else
        echo -e "  ${RED}✘ ${name}${NC} — ejecutar: ${cmd}"
    fi
done

fi  # TAX_ONLY

# ══════════════════════════════════════════════════════════════════
# RESUMEN
# ══════════════════════════════════════════════════════════════════

echo -e "\n${BLUE}── RESUMEN ───────────────────────────────────────────${NC}\n"
echo -e "  DB Root: ${DB_ROOT}"
echo ""

TOTAL_OK=0
TOTAL_MISSING=0

check_and_count() {
    local name="$1" dir="$2" marker="$3"
    if db_ok "$dir" "$marker"; then
        echo -e "  ${GREEN}✔${NC} $name"
        ((TOTAL_OK++))
    else
        echo -e "  ${RED}✘${NC} $name"
        ((TOTAL_MISSING++))
    fi
}

echo "  Fase TAX:"
check_and_count "Kraken2 PlusPF"    "${DB_ROOT}/k2_pluspf_20251015"     "hash.k2d"
check_and_count "Kaiju RefSeq"      "${DB_ROOT}/kaiju"                   "kaiju_db_refseq_ref.fmi"
check_and_count "Sylph GTDB"        "${DB_ROOT}/sylph"                   "gtdb-r226-c200-dbv1.syldb"
check_and_count "Sylph Fungi"       "${DB_ROOT}/sylph"                   "fungi-refseq-2024-07-25-c200-v0.3.syldb"
check_and_count "Sylph Viral"       "${DB_ROOT}/sylph"                   "imgvr_c200_v0.3.0.syldb"
check_and_count "ResFinder"         "${DB_ROOT}/kma_resfinder/resfinder_db"  "all.fsa"
check_and_count "FastQ Screen"      "${DB_ROOT}/FastqScreen"             "fastq_screen.conf"

if ! $TAX_ONLY; then
echo "  Fase MAG:"
check_and_count "GTDB-Tk r220"      "${DB_ROOT}/gtdbtk"      "release220/metadata/metadata.txt"
check_and_count "CheckM2"           "${DB_ROOT}/checkm2"     "CheckM2_database/uniref100.KO.1.dmnd"
check_and_count "Bakta"             "${DB_ROOT}/bakta"       "db"
check_and_count "eggNOG"            "${DB_ROOT}/eggnog"      "eggnog.db"
check_and_count "BUSCO"             "${DB_ROOT}/busco"       "lineages"
check_and_count "GUNC"              "${DB_ROOT}/gunc"        "gunc_db_progenomes2.1.dmnd"
check_and_count "geNomad"           "${DB_ROOT}/genomad"     "genomad_db"
check_and_count "CheckV"            "${DB_ROOT}/checkv"      "checkv-db-v1.5/README.txt"
check_and_count "AMRFinderPlus"     "${DB_ROOT}/amrfinderplus" "latest"
fi

echo ""
echo -e "  ${GREEN}${TOTAL_OK} OK${NC} | ${RED}${TOTAL_MISSING} faltan${NC}"

if [[ $TOTAL_MISSING -gt 0 ]]; then
    echo ""
    echo -e "  ${YELLOW}Para usar con el pipeline:${NC}"
    echo -e "    nextflow run main.nf --input data/... --db_root ${DB_ROOT} -profile gru"
else
    echo ""
    echo -e "  ${GREEN}¡Todas las bases de datos están listas!${NC}"
    echo -e "    nextflow run main.nf --input data/... --db_root ${DB_ROOT} -profile gru"
fi
echo ""
