#!/usr/bin/env bash

source /opt/conda/etc/profile.d/conda.sh

if [[ $# -ne 1 ]]; then
    echo "Usage: run_short.sh <mandatory_vars.txt>"
    exit 1
fi

CONF="$1"

if [[ ! -f "$CONF" ]]; then
    echo "[ERROR] Config file not found: $CONF"
    exit 1
fi

echo "[INFO] Using config file: $CONF"

# Load conda and activate environment
RAAVIOLI_SHORT_DIR="/opt/RAAVioli/RAAVioli_short"
cd "$RAAVIOLI_SHORT_DIR"

RAAVIOLI_SHORT="/opt/RAAVioli/RAAVioli_short/RAAVioli_short.sh"

if [[ ! -f "$RAAVIOLI_SHORT" ]]; then
    echo "[ERROR] RAAVioliShort.sh not found at $RAAVIOLI_SHORT"
    exit 1
fi

echo "[INFO] Running RAAVioli Short..."
bash "$RAAVIOLI_SHORT" "$CONF"

echo "[INFO] RAAVioli Short completed successfully."
