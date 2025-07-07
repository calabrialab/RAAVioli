for TAG in ${ASSOBCLIST[@]}; do
	##### ========================== REMOVING SEQUENCES THAT DON'T START WHITH FUSION  ================================= #####
	echo "REMOVING SEQUENCES THAT DON'T START WHITH FUSION"
	flexbar -r ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_1.fastq.gz -b ${FUSION_PRIMERS}  -bu -bt LTAIL -t ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_1_fusion

	awk 'NR%4==1 {print substr($1,2)}' ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_1_fusion_barcode_unassigned.fastq > ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_1_fusion.exclude.list

	zcat ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_1.fastq.gz  | \
	python3 ${RAAVIOLIDIR}/utils_scripts/fqextract.pureheader.v2.py ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_1_fusion.exclude.list exclude | pigz -f -c > ${TMPDIR}/bcmuxall/r1.no12.${TAG}.fastq.gz
	zcat ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_2.fastq.gz  | \
	python3 ${RAAVIOLIDIR}/utils_scripts/fqextract.pureheader.v2.py ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_1_fusion.exclude.list exclude | pigz -f -c > ${TMPDIR}/bcmuxall/r2.no12.${TAG}.fastq.gz

	#			rm ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_1_fusion*

	# # ##### =========================================================================================================== #####


	##### ================ Vector quantification TAG based - start ======================== #####
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference VECTOR genome"
	echo -e "R1: ${TMPDIR}/bcmuxall/r1.no12.${TAG}.fastq.gz \nR2: ${TMPDIR}/bcmuxall/r2.no12.${TAG}.fastq.gz"
	## --- paired-end reads ---
	# low seed value due to small genome
	bwa mem -k 14 -r 1 -v 1 -T ${BWA_MIN_ALN_LEN} -c 1 -R "@RG\tID:${b}_${k}\tSM:Vector\tCN:TIGET" -t ${MAXTHREADS} ${VECTORGENOME} <(zcat ${TMPDIR}/bcmuxall/r1.no12.${TAG}.fastq.gz ) | samtools view -F 260 -q ${mapQvec} -uS - | samtools sort - ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${TAG}.noLTRLC.sorted.md
	samtools index ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${TAG}.noLTRLC.sorted.md.bam ;
	# extract reads ID of the LTR mapping reads
	samtools view ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${TAG}.noLTRLC.sorted.md.bam | cut -f1 | sort | uniq > ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${TAG}.sorted.md.vector.list
	# fastqc -o ${OUTDIR_POOL_QUAL} --contaminants ${SEQCONTAM} -t ${MAXTHREADS} -f bam ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${TAG}.noLTRLC.merge.bam
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Extract no LV reads from raw data"
	### ${TMPDIR}/bcmuxall/r1.no12.${TAG}.noLTRLC.fastq.gz
	zcat ${TMPDIR}/bcmuxall/r1.no12.${TAG}.fastq.gz | python3 ${RAAVIOLIDIR}/utils_scripts/fqextract_pureheader.v3.py ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${TAG}.sorted.md.vector.list | pigz -f -c > ${TMPDIR}/bcmuxall/r1.no12.${TAG}.Vector.fastq.gz
	zcat ${TMPDIR}/bcmuxall/r2.no12.${TAG}.fastq.gz | python3 ${RAAVIOLIDIR}/utils_scripts/fqextract_pureheader.v3.py ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${TAG}.sorted.md.vector.list | pigz -f -c > ${TMPDIR}/bcmuxall/r2.no12.${TAG}.Vector.fastq.gz
	# remove
	# rm ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${TAG}.noLTRLC.sorted.md.bam
	pigz -f ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${TAG}.sorted.md.vector.list
done

#### ======================================= ALIGNMENT PAIRED ENDS ========= =======================================

for TAG in ${ASSOBCLIST[@]}; do
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference genome in PE and SE using BWA MEM"
	# test with  -P -M -c 2 -T 15 -> very stringent... too few reads in output
	# explaination: -c 1 -> keep only reads with 1 MEM; -P -> try to rescue some pairs with SW:: ATTENTION!!!!! This option will destroy all paired-end reads!!!!! Do NOT USE IT!; -T 15 -> min alignment score; -M -> flag short alignments (for Picard);
	bwa mem -k 18 -r 1 -M -v 1 -T 15 -c 1 -R "@RG\tID:${TAG}\tSM:${TAG}\tCN:TIGET.${DISEASE}.${PATIENT}.${POOL}" -t ${MAXTHREADS} ${GENOME} <(zcat ${TMPDIR}/bcmuxall/r1.no12.${TAG}.Vector.fastq.gz ) <(zcat ${TMPDIR}/bcmuxall/r2.no12.${TAG}.Vector.fastq.gz ) > ${TMPDIR}/sam/${TAG}.sam

	# create BAM and sort them
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating BAM and indexes (filter from here the dataset using only valid reads: mapped and primary)"
	samtools view -F 4 -uS ${TMPDIR}/sam/${TAG}.sam | samtools sort - ${TMPDIR}/bam/${TAG}.F4.sorted.md; # THIS CAN BE IMPROVED BY SAVING ONLY PRIMARY ALIGNMENT ON R2 WITH ITS SUPPLEMENTARY OR OTHERWISE BY FILTERING OUT THE READS NAME FROM THE R1 FINAL BAM
	samtools index ${TMPDIR}/bam/${TAG}.F4.sorted.md.bam
	samtools view -F 2308 -uS ${TMPDIR}/sam/${TAG}.sam | samtools sort - ${TMPDIR}/bam/${TAG}.sorted.md;
	# Removing Reads in chrM and chrUn
	samtools view -h ${TMPDIR}/bam/${TAG}.sorted.md.bam | awk '{if($3 != "chrM" && $3 != "chrUn"){print $0}}' | samtools view -Sb - > ${TMPDIR}/bam/${TAG}.sorted.md.cleaned.bam
	samtools index ${TMPDIR}/bam/${TAG}.sorted.md.cleaned.bam
	rm ${TMPDIR}/sam/${TAG}.sam ;
done


#### ======================================= FILTERING ================================================
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering data Branching Paired/Single-end mapped reads"
for TAG in ${ASSOBCLIST[@]}; do
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering data (Bamtools)"
	bamtools filter -in ${TMPDIR}/bam/${TAG}.sorted.md.cleaned.bam -isMapped true -isMateMapped true -isPaired true -isPrimaryAlignment true -mapQuality ">=${minmapQ}" -out ${TMPDIR}/bam/${TAG}.sorted.md.rel.pg.iss.bam
	samtools index ${TMPDIR}/bam/${TAG}.sorted.md.rel.pg.iss.bam
done
