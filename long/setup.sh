#!/bin/bash
DIR=$(pwd)
mkdir $DIR/bin
mkdir $DIR/logs
cd $DIR/bin
echo "Installing bwa"
git clone https://github.com/lh3/bwa.git
cd bwa
make > bwa.install.log 2>bwa.error.install.log
BWA=${DIR}/bin/bwa/bwa
mv bwa.*install.log $DIR/logs
cd ${DIR}

cd $DIR/bin
echo "Installing samtools"
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9
make > samtools.install.log 2>samtools.error.install.log
SAMTOOLS=${DIR}/bin/samtools-1.9/samtools
mv samtools.*install.log $DIR/logs
cd ${DIR}


cd $DIR/bin
echo "Installing bamtools"
git clone https://github.com/pezmaster31/bamtools.git
cd bamtools
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${DIR}/bin/bamtools .. > bamtools.install.log 2>bamtools.error.install.log
make >> bamtools.install.log 2>>bamtools.error.install.log
make install >> bamtools.install.log 2>>bamtools.error.install.log
mv bamtools.*install.log $DIR/logs
BAMTOOLS=${DIR}/bin/bamtools/bin/bamtools


cd $DIR/bin
echo "Installing bedtools"
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make > bedtools.install.log 2>bedtools.error.install.log
BEDTOOLS=${DIR}/bin/bedtools2/bin/bedtools
mv bedtools.*install.log $DIR/logs
cd ${DIR}


cd $DIR/bin
echo "Installing fastq_to_fasta from fastx_toolkit"  
mkdir fastx_toolkit
cd fastx_toolkit
wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
tar -xf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
FASTQ_TO_FASTA=${DIR}/bin/fastx_toolkit/bin/fastq_to_fasta
cd ${DIR}





touch ${DIR}/config.txt
echo "BWA=${BWA}" >> config.txt
echo "SAMTOOLS=${SAMTOOLS}" >> config.txt
echo "BAMTOOLS=${BAMTOOLS}" >> config.txt
echo "BEDTOOLS=${BEDTOOLS}" >> config.txt
echo "FASTQ_TO_FASTA=${FASTQ_TO_FASTA}" >> config.txt
echo "VARIABLES_VIRAL=${DIR}/variables_viral" >> config.txt
echo "VARIABLES_MIXED=${DIR}/variables_mixed" >> config.txt
echo "VARIABLES_STEPR=${DIR}/variables_rscript" >> config.txt
echo "FQEXTRACT=${DIR}/scripts/fqextract.pureheader.py" >> config.txt
echo "FASTA_TO_CSV=${DIR}/scripts/fasta_to_csv.rb" >> config.txt
echo "Done"
