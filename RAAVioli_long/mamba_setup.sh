#!/bin/bash
set -euo pipefail

# --- SETTINGS ---
ENV_NAME="RAAVioliLong_env"
MAMBA_ENV="mambarav_env"        # lightweight env only for mamba (created if needed)
YML_FILE="raaviolilong_env.yml"
DIR="$(pwd)"
LOG_DIR="$DIR/logs"
CONFIG_FILE="$DIR/config.txt"

# --- PRECHECKS ---
mkdir -p "$LOG_DIR"

if ! command -v conda &> /dev/null; then
  echo "[ERROR] Conda not found. Please install Miniconda or Mambaforge first and re-run this script."
  exit 1
fi

# --- STEP 1: LOCATE OR INSTALL MAMBA ---
if command -v mamba &> /dev/null; then
  echo "[INFO] Using system mamba: $(which mamba)"
  MAMBA_CMD="mamba"
else
  echo "[WARN] Mamba not found. Creating helper environment '$MAMBA_ENV'..."
  if conda env list | grep -q "^$MAMBA_ENV\s"; then
    echo "[INFO] Helper environment '$MAMBA_ENV' already exists."
  else
    conda create -y -n "$MAMBA_ENV" -c conda-forge python=3.10 mamba
  fi

  MAMBA_CMD="conda run -n $MAMBA_ENV mamba"
  echo "[INFO] Using mamba from helper environment '$MAMBA_ENV'."
fi

# --- STEP 2: CREATE OR UPDATE MAIN ENVIRONMENT ---
if conda env list | grep -q "^$ENV_NAME\s"; then
  echo "[INFO] Environment '$ENV_NAME' already exists. Updating it..."
  if [[ -f "$YML_FILE" ]]; then
    $MAMBA_CMD env update -y -n "$ENV_NAME" -f "$YML_FILE"
  else
    echo "[WARN] Environment file '$YML_FILE' not found. Skipping YAML update."
  fi
else
  echo "[INFO] Creating new environment '$ENV_NAME'..."
  $MAMBA_CMD create -y -n "$ENV_NAME" -c conda-forge python=3.9
  if [[ -f "$YML_FILE" ]]; then
    $MAMBA_CMD env update -y -n "$ENV_NAME" -f "$YML_FILE"
  fi
fi

# --- STEP 3: ACTIVATE MAIN ENVIRONMENT ---
eval "$(conda shell.bash hook)"
set +u
conda activate "$ENV_NAME"
set -u

# --- FIX LOCALE (for R LC_CTYPE errors) ---
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# --- STEP 4: INSTALL R PACKAGES VIA CONDA ---
echo "[INFO] Installing R base packages from conda-forge..."
$MAMBA_CMD install -y -n "$ENV_NAME" \
  -c conda-forge -c bioconda \
  r-base r-optparse r-sqldf r-matrix \
  bioconductor-delayedarray bioconductor-summarizedexperiment

# --- STEP 5: INSTALL CRAN & BIOC PACKAGES IN R ---
echo "[INFO] Installing CRAN and Bioconductor packages inside R..."
Rscript - <<'EOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))

packages <- c("optparse", "tools", "sqldf")

to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

if (!requireNamespace("GenomicAlignments", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("GenomicAlignments", update = FALSE, ask = FALSE)
}
EOF

# --- STEP 6: LOG TOOL VERSIONS ---
{
  echo "[SETUP LOG - $(date)]"
  echo "BWA: $(bwa 2>&1 | head -n1 || echo 'not found')"
  echo "SAMTOOLS: $(samtools --version | head -n1 || echo 'not found')"
  echo "BAMTOOLS: $(bamtools --version 2>&1 | head -n1 || echo 'not found')"
  echo "BEDTOOLS: $(bedtools --version 2>/dev/null || echo 'not found')"
  echo "R: $(R --version | head -n1 || echo 'not found')"
  echo "Python: $(python --version 2>&1)"
} > "$LOG_DIR/versions.log"

# --- STEP 7: MAKE SCRIPTS EXECUTABLE ---
chmod +x "$DIR/scripts/fqextract.pureheader.v3.py" 2>/dev/null || true
chmod +x "$DIR/scripts/fasta_to_csv.rb" 2>/dev/null || true

# --- STEP 8: GENERATE CONFIG FILE ---
{
  echo "BWA=bwa"
  echo "SAMTOOLS=samtools"
  echo "BAMTOOLS=bamtools"
  echo "BEDTOOLS=bedtools"
  echo "FASTQ_TO_FASTA=$DIR/scripts/fastq_to_fasta.tiget.v3.py"
  echo "FQEXTRACT=$DIR/scripts/fqextract.pureheader.v3.py"
  echo "FASTA_TO_CSV=$DIR/scripts/fasta_to_csv.rb"
} > "$CONFIG_FILE"

# --- SUCCESS MESSAGE ---
echo ""
echo "[SUCCESS] RAAVioliLongR environment setup complete."
echo "  - Environment name: $ENV_NAME"
if [[ "$MAMBA_CMD" != "mamba" ]]; then
  echo "  - Helper mamba env: $MAMBA_ENV"
fi
echo "  - Config file: $CONFIG_FILE"
echo "  - Log file: $LOG_DIR/versions.log"
echo ""
echo "To activate your environment, run:"
echo "    conda activate $ENV_NAME"
