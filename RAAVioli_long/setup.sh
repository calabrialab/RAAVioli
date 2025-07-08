#!/bin/bash
set -euo pipefail

# Set working directory
DIR="$(pwd)"
mkdir -p "$DIR/logs"

# Check if conda is installed
if ! command -v conda &> /dev/null; then
  echo "Conda not found. Please install Miniconda or Anaconda first and re-run this script."
  exit 1
fi

# Conda environment name
ENV_NAME="RAAVioliLong_env"

echo "Creating conda environment '$ENV_NAME' "

# Create environment with specific channels and versions
conda create -y -n "$ENV_NAME" \
  -c conda-forge \
  -c bioconda \
  -c defaults \
  bwa=0.7.17 \
  samtools=1.9 \
  bamtools=2.5.1 \
  bedtools=2.29.1

# Activate environment
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

# Log tool versions
{
  echo "BWA: $(bwa 2>&1 | head -n1)"
  echo "SAMTOOLS: $(samtools --version | head -n1)"
  echo "BAMTOOLS: $(bamtools --version 2>&1 | head -n1)"
  echo "BEDTOOLS: $(bedtools --version)"
} > "$DIR/logs/versions.log"

chmod +x "$DIR/scripts/fqextract.pureheader.v3.py"
chmod +x "$DIR/scripts/fasta_to_csv.rb"

# Generate config file
CONFIG="$DIR/config.txt"
{
  echo "BWA=bwa"
  echo "SAMTOOLS=samtools"
  echo "BAMTOOLS=bamtools"
  echo "BEDTOOLS=bedtools"
  echo "FASTQ_TO_FASTA=$DIR/scripts/fastq_to_fasta.tiget.v3.py"
  echo "FQEXTRACT=$DIR/scripts/fqextract.pureheader.v3.py"
  echo "FASTA_TO_CSV=$DIR/scripts/fasta_to_csv.rb"
} > "$CONFIG"

echo "Environment and config file setup complete."
echo "To activate your environment: conda activate $ENV_NAME"
