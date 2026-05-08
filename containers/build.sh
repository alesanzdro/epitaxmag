#!/bin/bash
# Build script for the EpiTaxMAG Tools image.
#
# OPTION 1: Build with Singularity (direct, no Docker required)
#   bash containers/build.sh singularity
#
# OPTION 2: Build with Docker and convert to Singularity
#   bash containers/build.sh docker
#
# The resulting image is written to containers/epitaxmag-tools-1.0.0.sif

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
        echo "Usage: bash containers/build.sh [singularity|docker]"
        exit 1
        ;;
esac

echo ""
echo "Image created: ${SIF_FILE}"
echo ""
echo "Verify:"
echo "  singularity exec ${SIF_FILE} python3 -c 'import pandas, plotly; print(\"OK\")'"
echo "  singularity exec ${SIF_FILE} minimap2 --version"
echo "  singularity exec ${SIF_FILE} kma -v"
