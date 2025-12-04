#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: run_long.sh <RAAVioli_long.sh arguments>"
    echo ""
    echo "Example:"
    echo "  run_long.sh -i /opt/config/sample_label.tsv -t 8 -V /opt/genomes/vector.fa -M /opt/genomes/mixed.fa -a /opt/annot/custom.gtf -o /opt/output -c vars_mixed -w vars_viral -y vars_rscript"
    exit 1
fi


RAAVIOLI_LONG="/opt/RAAVioli/RAAVioli_long/RAAVioli_long.sh"

if [[ ! -f "$RAAVIOLI_LONG" ]]; then
    echo "[ERROR] RAAVioli_long.sh not found at $RAAVIOLI_LONG"
    exit 1
fi

echo "[INFO] Running RAAVioli Long with arguments:"
echo "      $@"

bash "$RAAVIOLI_LONG" "$@"

echo "[INFO] RAAVioli Long completed successfully."
