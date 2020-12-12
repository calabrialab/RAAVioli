#!/bin/bash

echo "Checking if dependencies are already installed"
DIR=$(pwd)
mkdir $DIR/bin

BWA=$(which bwa)
if [ -z "$BWA" ]
then
    cd $DIR/bin
    echo "Installing bwa"
    git clone https://github.com/lh3/bwa.git
    cd bwa
    make
    BWA=${DIR}/bin/bwa/bwa
    cd ${DIR}
fi

SAMTOOLS=$(which samtools)
if [ -z "$SAMTOOLS" ]
then
    cd $DIR/bin
    echo "Installing samtools"
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    tar -vxjf samtools-1.9.tar.bz2
    cd samtools-1.9
    make
    BWA=${DIR}/bin/samtools-1.9/samtools
    cd ${DIR}
fi

BAMTOOLS=$(which bamtools)
if [ -z "$BAMTOOLS" ]
then
    echo "Installing bamtools"
    sudo apt-get install bamtools
    BAMTOOLS=$(which bamtools)
fi

BEDTOOLS=$(which bedtools)
if [ -z "$BEDTOOLS" ]
then
    cd $DIR/bin
    echo "Installing bedtools"
    wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
    tar -zxvf bedtools-2.29.1.tar.gz
    cd bedtools2
    make
    BEDTOOLS=${DIR}/bin/bedtools2/bin/bedtools
    cd ${DIR}
fi
    
FASTQ_TO_FASTA=$(which fastq_to_fasta)
if [ -z "$FASTQ_TO_FASTA" ]  
then
    cd $DIR/bin
    echo "Installing fastq_to_fasta from fastx_toolkit"  
    mkdir fastx_toolkit
    cd fastx_toolkit
    wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
    tar -xf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
    FASTQ_TO_FASTA=${DIR}/bin/fastx_toolkit/bin/fastq_to_fasta
    cd ${DIR}
fi

BAMCOVERAGE=$(which bamCoverage)
if [ -z "$BAMCOVERAGE" ]  
then
    echo "Installing bamCoverage froom deepTools"  
    pip install deepTools --user
    BAMCOVERAGE=$(which bamCoverage)
fi

touch ${DIR}/config.txt
echo "BWA=${BWA}" >> config.txt
echo "SAMTOOLS=${SAMTOOLS}" >> config.txt
echo "BAMTOOLS=${BAMTOOLS}" >> config.txt
echo "BEDTOOLS=${BEDTOOLS}" >> config.txt
echo "FASTQ_TO_FASTA=${FASTQ_TO_FASTA}" >> config.txt
echo "BAMCOVERAGE=${BAMCOVERAGE}" >> config.txt
echo "VARIABLES_VIRAL=${DIR}/variables_viral" >> config.txt
echo "VARIABLES_REFERENCE=${DIR}/variables_reference" >> config.txt
echo "VARIABLES_MIXED=${DIR}/variables_mixed" >> config.txt
echo "FQEXTRACT=${DIR}/scripts/fqextract.pureheader.py" >> config.txt
echo "FASTA_TO_CSV=${DIR}/scripts/fasta_to_csv.rb" >> config.txt
echo "Done"
