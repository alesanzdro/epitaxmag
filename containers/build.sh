#!/bin/bash
# Build de la imagen EpiTaxMAG Tools
#
# OPCION 1: Build con Singularity (directo, sin Docker)
#   bash containers/build.sh singularity
#
# OPCION 2: Build con Docker + convertir a Singularity
#   bash containers/build.sh docker
#
# La imagen resultante se guarda en containers/epitaxmag-tools-1.0.0.sif

set -euo pipefail

VERSION="1.0.0"
IMAGE_NAME="epitaxmag-tools"
SIF_FILE="containers/${IMAGE_NAME}-${VERSION}.sif"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

MODE="${1:-singularity}"

echo ""
echo "  Building ${IMAGE_NAME}:${VERSION}"
echo "  Mode: ${MODE}"
echo ""

case "$MODE" in
    singularity)
        echo "Building with Singularity (--fakeroot)..."
        singularity build --fakeroot "$SIF_FILE" "${SCRIPT_DIR}/epitaxmag-tools.def"
        ;;
    docker)
        echo "Building Docker image..."
        docker build -t "${IMAGE_NAME}:${VERSION}" -f "${SCRIPT_DIR}/Dockerfile.epitaxmag-tools" .
        echo "Exporting Docker image to tar..."
        docker save "${IMAGE_NAME}:${VERSION}" -o "${SCRIPT_DIR}/${IMAGE_NAME}-${VERSION}.tar"
        echo "Converting tar to Singularity SIF..."
        singularity build "$SIF_FILE" "docker-archive://${SCRIPT_DIR}/${IMAGE_NAME}-${VERSION}.tar"
        rm -f "${SCRIPT_DIR}/${IMAGE_NAME}-${VERSION}.tar"
        ;;
    *)
        echo "Uso: bash containers/build.sh [singularity|docker]"
        exit 1
        ;;
esac

echo ""
echo "Imagen creada: ${SIF_FILE}"
echo ""
echo "Verificar:"
echo "  singularity exec ${SIF_FILE} python3 -c 'import pandas, plotly; print(\"OK\")'"
echo "  singularity exec ${SIF_FILE} minimap2 --version"
echo "  singularity exec ${SIF_FILE} kma -v"
