#!/usr/bin/env bash

# === FUNCTIONS === #
usage() {
    echo "Usage: $0 -g <host_genome.fa> -v <vector_genome.fa> -o <output_dir> -n <output_name.fa>"
    echo "Example: $0 -g /data/hg19.fa -v /data/aav.fa -o /data/mixed/ -n mixedgenome.fa"
    exit 1
}

timestamp() {
    date +'%Y-%m-%d %H:%M:%S'
}

log() {
    echo "[AP] ============ <$(timestamp)> [TIGET] $1 ============"
}

# === PARSE ARGUMENTS === #
while getopts ":g:v:o:n:" opt; do
    case ${opt} in
        g) HOSTGENOME="$OPTARG" ;;
        v) VECTORGENOME="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        n) OUTPUT_NAME="$OPTARG" ;;
        *) usage ;;
    esac
done

# === CHECK ARGUMENTS === #
if [[ -z "$HOSTGENOME" || -z "$VECTORGENOME" || -z "$OUTPUT_DIR" || -z "$OUTPUT_NAME" ]]; then
    usage
fi

# === CHECK INPUT FILES EXIST === #
if [[ ! -f "$HOSTGENOME" ]]; then
    echo "[ERROR] Host genome not found at: $HOSTGENOME"
    exit 1
fi

if [[ ! -f "$VECTORGENOME" ]]; then
    echo "[ERROR] Vector genome not found at: $VECTORGENOME"
    exit 1
fi

# === CHECK VECTOR SEQUENCE NAME === #
CHR_NAME=$(grep "^>" "$VECTORGENOME" | head -n 1 | sed 's/>//')
if [[ "$CHR_NAME" != "chrV" ]]; then
    echo "[ERROR] Vector genome must have sequence name 'chrV' (found: '$CHR_NAME')"
    exit 1
fi

# === CREATE OUTPUT DIR IF NEEDED === #
if [[ ! -d "$OUTPUT_DIR" ]]; then
    log "Output directory does not exist, creating: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR" || { echo "[ERROR] Failed to create output directory."; exit 1; }
fi

# === DEFINE OUTPUT FILE === #
MIXEDGENOME="${OUTPUT_DIR%/}/${OUTPUT_NAME}"

# === CHECK IF MIXED GENOME ALREADY EXISTS === #
if [[ -f "$MIXEDGENOME" ]]; then
    echo "[ERROR] Mixed genome already exists at: $MIXEDGENOME"
    echo "        Please delete it manually before re-running this script."
    exit 1
fi

# === CREATE MIXED GENOME === #
log "Creating mixed genome at ${MIXEDGENOME}"
cp "$HOSTGENOME" "$MIXEDGENOME"
cat "$VECTORGENOME" >> "$MIXEDGENOME"

# === INDEX MIXED GENOME === #
if ! command -v bwa &>/dev/null; then
    echo "[ERROR] bwa not found in PATH. Please install or load bwa before running."
    exit 1
fi

log "Indexing mixed genome"
bwa index -a bwtsw "$MIXEDGENOME"

log "Done"
