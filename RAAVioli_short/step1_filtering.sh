
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] New TMP files check"
ls -lh ${R1_FASTQ}
ls -lh ${R2_FASTQ}

##### ================ RAW DATA QUALITY ======================== #####
# 0. quality check ...! NOT WORKING FINE IIF FASTQ.GZ COMES FROM ILLUMINA MACHINE/SW...!!! YOU MUST CONVERT THEM ALL: gunzip first and then gzip them all again (their compression is not standard)
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Quality check -> FastQC reports and FastX stats"

quack -1 ${R1_FASTQ} -2 ${R2_FASTQ} -n ${DISEASE}-${PATIENT}-${POOL} > ${OUTDIR_POOL_QUAL}/quack-analysis.svg

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Quality check -> FastQ-Screen reports across genomes"
fastq_screen ${R1_FASTQ} ${R2_FASTQ} --outdir ${OUTDIR_POOL_QUAL} --subset 500000 --aligner bwa --force

# FASTA of RANDOM TAGs
if ! [ ${FASTQ_QF} = "no_quality_filter" ]; then
	if [ -z "$Q_FILTER" ]; then
        fastp -i ${R1_FASTQ} -I ${R2_FASTQ} -o ${TMPDIR}/out.R1.fq.gz -O ${TMPDIR}/out.R2.fq.gz --disable_quality_filtering --fix_mgi_id -w ${MAXTHREADS} --disable_adapter_trimming
		## remove low quality R1-R2 reads with fastq_qf program
		fastq_qf -a ${R1_FASTQ} -b ${R2_FASTQ} -o ${TMPDIR} -t ${MAXTHREADS} -m ${FASTQ_QF}
		R1_FASTQ="${TMPDIR}/QF.${R1_NAME}";
		R2_FASTQ="${TMPDIR}/QF.${R2_NAME}";
		QUALITY_PASSED=$((`zcat ${R1_FASTQ} | wc -l `/4));
		rm ${TMPDIR}/${R1_NAME} ${TMPDIR}/${R2_NAME};
	else
        fastp -i ${R1_FASTQ} -I ${R2_FASTQ} -o ${TMPDIR}/out.R1.fq.gz -O ${TMPDIR}/out.R2.fq.gz --disable_quality_filtering --fix_mgi_id -w ${MAXTHREADS} --disable_adapter_trimming
		## remove low quality R1-R2 reads with fastq_qf program
		fastq_qf -a ${R1_FASTQ} -b ${R2_FASTQ} -o ${TMPDIR} -t ${MAXTHREADS} -m ${FASTQ_QF} -f ${Q_FILTER}
		R1_FASTQ="${TMPDIR}/QF.${R1_NAME}";
		R2_FASTQ="${TMPDIR}/QF.${R2_NAME}";
		QUALITY_PASSED=$((`zcat ${R1_FASTQ} | wc -l `/4));
		rm ${TMPDIR}/${R1_NAME} ${TMPDIR}/${R2_NAME};
	fi
fi

if [ ${FASTQ_QF} == "no_quality_filter" ]; then
	fastp -i ${R1_FASTQ} -I ${R2_FASTQ} -o ${TMPDIR}/out.R1.fq.gz -O ${TMPDIR}/out.R2.fq.gz --average_qual 31 --fix_mgi_id -w ${MAXTHREADS} --disable_adapter_trimming
	QUALITY_PASSED=$((`zcat ${TMPDIR}/out.R1.fq.gz | wc -l `/4));
	R1_FASTQ="${TMPDIR}/out.R1.fq.gz";
	R2_FASTQ="${TMPDIR}/out.R2.fq.gz";
	mv fastp.html fastp.json ${OUTDIR_POOL_QUAL}
fi

# ##### ================ PHIX quantification START ======================== #####
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] PHIX Alignment to reference genome each single pair"

### NEW version (without sam output)
#bwa mem -k 16 -r 1 -M -v 1 -T 15 -t ${MAXTHREADS} ${PHIXGENOME} <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} ) | samtools view -F 2308 -q 25 -f 35 -uS - > ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.phix.PE.bam;

bwa mem -k 16 -r 1 -M -v 1 -T 15 -t "${MAXTHREADS}" "${PHIXGENOME}" <(zcat "${R1_FASTQ}") <(zcat "${R2_FASTQ}") \
| samtools view -F 2308 -f 35 -q 25 -O BAM,uncompressed -o "${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.phix.PE.bam"


echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering raw data"
samtools view ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.phix.PE.bam | cut -f 1 > ${TMPDIR}/PHIX.header.list;
sort --parallel=5 ${TMPDIR}/PHIX.header.list > ${TMPDIR}/PHIX.header.sorted.list;
rm ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.phix.PE.bam ${TMPDIR}/PHIX.header.list;
## rename files
BNAME_R1=`basename ${R1_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
BNAME_R2=`basename ${R2_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
zcat ${R1_FASTQ} | fqreverseextract_pureheader ${TMPDIR}/PHIX.header.sorted.list | pigz -f -c > ${TMPDIR}/${BNAME_R1}_R1_nophix.fastq.gz &
zcat ${R2_FASTQ} | fqreverseextract_pureheader ${TMPDIR}/PHIX.header.sorted.list | pigz -f -c > ${TMPDIR}/${BNAME_R2}_R2_nophix.fastq.gz &
wait
# checking VARS
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Checking Variables ";
PHIX_MAPPING=`wc -l ${TMPDIR}/PHIX.header.sorted.list | cut -d' ' -f1 ` ;
##PHIX_MAPPING="99999999" ;
## give new gzipped files without PHIX, change variables
rm ${R1_FASTQ} ${R2_FASTQ} ${TMPDIR}/PHIX.header.sorted.list
R1_FASTQ="${TMPDIR}/${BNAME_R1}_R1_nophix.fastq.gz" # BNAME_R1=`basename ${R1_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
R2_FASTQ="${TMPDIR}/${BNAME_R2}_R2_nophix.fastq.gz"
# ##### ================ PHIX quantification END  ======================== #####


##### ================ PLASMIDS quantification POOL based - start ======================== #####
BNAME_R1=`basename ${R1_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
BNAME_R2=`basename ${R2_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;

ATLEASTAPLASMID=0;
for PLASMID_FILE in ${PLASMID_DIR}/plasmids.*.fa; do
    if [[ "$PLASMID_FILE" != "${PLASMID_DIR}/plasmids.*.fa" ]]; then
		PLASMID=`basename ${PLASMID_FILE} | sed 's/.fa//g' | cut -d'.' -f2`
		echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference ${PLASMID} genome"
		#bwa mem -k 14 -r 1 -v 1 -T 15 -c 1 -t ${MAXTHREADS} ${PLASMID_FILE} <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} ) | samtools view -F 2308 -q 20 -uS - | samtools sort - ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${PLASMID}.sorted
		bwa mem -k 14 -r 1 -v 1 -T 15 -c 1 -t "${MAXTHREADS}" "${PLASMID_FILE}" <(zcat "${R1_FASTQ}") <(zcat "${R2_FASTQ}") \
    | samtools view -F 2308 -q 20 -O BAM,uncompressed \
    | samtools sort -O BAM,uncompressed -o "${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${PLASMID}.sorted.bam"

		samtools view ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${PLASMID}.sorted.bam | cut -f1 > ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${PLASMID}.list
		ATLEASTAPLASMID=1;
    fi
done

if [[ "$ATLEASTAPLASMID" -eq 1 ]]; then
	sort --parallel=5 -u ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.*.list > ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.plasmids.list

	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Extract no PLASMIDS reads from raw data"
	### ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.noLTRLC.fastq.gz
	zcat ${R1_FASTQ} | fqreverseextract_pureheader ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.plasmids.list | pigz -f -c > ${TMPDIR}/${BNAME_R1}_noPlasmids.fastq.gz &
	zcat ${R2_FASTQ} | fqreverseextract_pureheader ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.plasmids.list | pigz -f -c > ${TMPDIR}/${BNAME_R2}_noPlasmids.fastq.gz &
	wait
	# var renaming
	rm ${R1_FASTQ} ${R2_FASTQ}
	R1_FASTQ="${TMPDIR}/${BNAME_R1}_noPlasmids.fastq.gz"
	R2_FASTQ="${TMPDIR}/${BNAME_R2}_noPlasmids.fastq.gz"
	# file compression
	pigz -f ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.plasmids.list
	pigz -f ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.*.list
fi
# get the number of remaining sequences
RAW_NO_PLASMID=$((`zcat ${R1_FASTQ} | wc -l | cut -d' ' -f1 `/4)) ;
##### ================ PLASMIDS quantification POOL based - end ======================== #####
