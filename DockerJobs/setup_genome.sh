#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo ""
    echo "Usage: setup_genome.sh --host <host.fa> --vector <vector.fa> --outdir <output_dir> [--name <mixed_name.fa>]"
    echo ""
    echo "Example:"
    echo "  setup_genome.sh --host /opt/host/hg19.fa --vector /opt/vector/vector.fa --outdir /opt/genomes/mixed --name mixed.fa"
    echo ""
    exit 1
}

HOST=""
VECTOR=""
OUTDIR=""
NAME="mixed.fa"   # default

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        --host) HOST="$2"; shift 2;;
        --vector) VECTOR="$2"; shift 2;;
        --outdir) OUTDIR="$2"; shift 2;;
        --name) NAME="$2"; shift 2;;
        *) usage;;
    esac
done

# Validate
if [[ -z "$HOST" || -z "$VECTOR" || -z "$OUTDIR" ]]; then
    echo "[ERROR] Missing required arguments."
    usage
fi

if [[ ! -s "$HOST" ]]; then
    echo "[ERROR] Host genome file not found or empty: $HOST"
    exit 1
fi

if [[ ! -s "$VECTOR" ]]; then
    echo "[ERROR] Vector genome file not found or empty: $VECTOR"
    exit 1
fi

mkdir -p "$OUTDIR"

MIXER="/opt/RAAVioli/RAAVioli_short/utils_scripts/create_mixed_genome.sh"

echo "[INFO] Creating mixed genome using RAAVioli script:"
echo "  Host:   $HOST"
echo "  Vector: $VECTOR"
echo "  Outdir: $OUTDIR"
echo "  Output name: $NAME"
echo ""

bash "$MIXER" \
    -g "$HOST" \
    -v "$VECTOR" \
    -o "$OUTDIR" \
    -n "$NAME"

echo ""
echo "[INFO] Mixed genome created successfully:"
echo "  $OUTDIR/$NAME"
