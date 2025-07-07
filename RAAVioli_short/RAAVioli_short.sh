RAAVIOLIDIR=$(pwd)
echo $RAAVIOLIDIR
# Check if the file is provided as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <mandatoryvars.txt>"
    exit 1
fi

mandatoryvars=$1

source $mandatoryvars

OUTDIR_MERGE_BAM="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/bam";
OUTDIR_POOL_BAM="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/bam/"${POOL};
OUTDIR_MERGE_BED="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/bed";
OUTDIR_POOL_BED="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/bed/"${POOL};
OUTDIR_MERGE_QUAL="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/quality";
OUTDIR_POOL_QUAL="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/quality/"${POOL};
OUTDIR_MERGE_QUANTIF="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/quantification";
OUTDIR_POOL_QUANTIF="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/quantification/"${POOL};
OUTDIR_MERGE_ISS="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/iss";
OUTDIR_POOL_ISS="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/iss/"${POOL};
OUTDIR_POOL_ISS_REPEATS="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/iss/"${POOL}/repeats;
OUTDIR_MERGE_BCMUXALL="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/bcmuxall";
OUTDIR_POOL_BCMUXALL="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/bcmuxall/"${POOL};
OUTDIR_MERGE_MATRIX="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/matrix/";
OUTDIR_POOL_MATRIX="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/matrix/"${POOL};


mkdir ${NGSWORKINGPATH}
mkdir ${NGSWORKINGPATH}/${DISEASE}
mkdir ${NGSWORKINGPATH}/${DISEASE}/${PATIENT}
mkdir ${OUTDIR_MERGE_MATRIX};
mkdir ${OUTDIR_POOL_MATRIX};
mkdir ${OUTDIR_MERGE_BAM};
mkdir ${OUTDIR_MERGE_BED};
mkdir ${OUTDIR_MERGE_QUAL};
mkdir ${OUTDIR_POOL_QUAL};
mkdir ${OUTDIR_MERGE_QUANTIF};
mkdir ${OUTDIR_POOL_QUANTIF};
mkdir ${OUTDIR_MERGE_ISS};
mkdir ${OUTDIR_POOL_ISS};
mkdir ${OUTDIR_POOL_ISS_REPEATS};
mkdir ${OUTDIR_MERGE_BCMUXALL};
mkdir ${OUTDIR_POOL_BCMUXALL};

mkdir ${OUTDIR_POOL_BAM};
mkdir ${OUTDIR_POOL_BED};
mkdir ${TMPDIR} ;

mkdir ${TMPDIR}/bcmuxall/
mkdir ${TMPDIR}/bed/
mkdir ${TMPDIR}/bam/
mkdir ${TMPDIR}/sam/
mkdir ${TMPDIR}/pools/
mkdir ${TMPDIR}/iss/
mkdir ${TMPDIR}/matrix/

mkdir ${TMPDIR}/stats/


OUTDIR_VECTOR=${TMPDIR}/pools
OUTDIR_VECTOR_POOL=${TMPDIR}/pools


# checking vars
RUN_ID=`date +"%Y%m%d%H%M%S"`
RUN_NAME="${DISEASE}|${PATIENT}|${POOL}"

BCLTRf=`cut "${BARCODE_LTR}" -f1`;
BCLTR=(${BCLTRf}) ;
BCLCf=`cut "${BARCODE_LC}" -f1`;
BCLC=(${BCLCf}) ;



# Specify the target string to search
target_string="TagID"

# Use awk to find the index of the target string in the header
target_index=$(awk -F'\t' -v target="$target_string" 'NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i == target) {
            print i
            exit
        }
    }
    print -1
}' ${ASSOCIATIONFILE})

if [ $target_index -ne -1 ]; then
    echo "Index of '$target_string' in the header: $target_index"
else
    echo "Error: '$target_string' not found in the header"
    exit 1
fi

ASSOBCLIST=(`cut -f${target_index} ${ASSOCIATIONFILE} | tail -n+2 | sort | uniq `); ## barcode list (as first colun of the association fTIGETile!) user defined!



##### ================ COPY DATA INTO TMP DIR ================== #####
# IF step1 performed
if [ -n "$step1_vars_file" ]; then
	source $step1_vars_file
	echo $R1_FASTQ
	cp ${R1_FASTQ} ${TMPDIR}/
	cp ${R2_FASTQ} ${TMPDIR}/

	# identify base names
	R1_NAME="`basename ${R1_FASTQ}`";
	R2_NAME="`basename ${R2_FASTQ}`";
	# now reassign raw input vars
	echo $R1_NAME
	R1_FASTQ=${TMPDIR}/${R1_NAME};
	R2_FASTQ=${TMPDIR}/${R2_NAME};
	RAW_READS=$((`zcat ${R1_FASTQ} | wc -l `/4));
	source step1_filtering.sh
fi

#HERE BEFORE WE HAVE TO CHECK THAT IF STEP1 HAS BEEN PERFORMED IT MUST PERFORM STEP2 OTHERWISE STEP3 will not work
if [ -n "$step2_vars_file" ]; then
	source $step2_vars_file
	source step2_demux.sh
fi

if [ -n "$step3_vars_file" ]; then
  # IF step1 and 2 not performed
	if [ -z "$step2_vars_file" ] && [ -z "$step1_vars_file" ]; then
		for TAG in ${ASSOBCLIST[@]}; do
		##### ==========================  MOVING READS IN THE BCMUXALL DIR IF STEP 1 AND 2 NOT DONE ================================= #####
			echo "MOVING READS IN THE BCMUXALL DIR"
			cp ${R1_FASTQ}/${TAG}.r1.fastq.gz  ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_1.fastq.gz
			cp ${R2_FASTQ}/${TAG}.r2.fastq.gz  ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_2.fastq.gz
    	done
	fi
	source $step3_vars_file
	source step3_alignment.sh
fi



if [ -n "$step4_vars_file" ]; then
  source $step4_vars_file
  source step4_is_identification.sh
fi




for TAG in ${ASSOBCLIST[@]}; do
	# checking VARS
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Checking Variables ";
	LV_MAPPING_RESULTS=`samtools flagstat ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${TAG}.noLTRLC.sorted.md.bam ` ; # ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${TAG}.noLTRLC.sorted.md.bam
	LV_MAPPING_PP=$((`echo ${LV_MAPPING_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
	LV_MAPPING_ST=$((`echo ${LV_MAPPING_RESULTS} | grep ' singletons ' | cut -d' ' -f48`)) ;
	LV_MAPPING_OVERALL=$((${LV_MAPPING_ST}+${LV_MAPPING_PP})) ;
	LV_MAPPING_MAPPED=$((`echo ${LV_MAPPING_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
	LV_MAPPED=`zcat ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${TAG}.sorted.md.vector.list.gz | wc -l | cut -d' ' -f1 `

	PLASMID_MAPPED_BYPOOL=`zcat ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.plasmids.list.gz | wc -l | cut -d' ' -f1 `

	BARCODE_MUX=$((`zcat ${TMPDIR}/bcmuxall/r1.no12.${TAG}.fastq.gz | wc -l `/4)) ;


	BWA_RESULTS=`samtools flagstat ${TMPDIR}/bam/${TAG}.sorted.md.cleaned.bam `
	BWA_INPUT_PAIRED=$((`zcat ${TMPDIR}/bcmuxall/r1.no12.${TAG}.Vector.fastq.gz  | wc -l | cut -d' ' -f1 `/4)) ;

	BWA_INPUT=$((${BWA_INPUT_PAIRED})) ;
	BWA_MAPPED=$((`echo ${BWA_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
	BWA_MAPPED_PP=$((`echo ${BWA_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
	BWA_MAPPED_ST=$((`echo ${BWA_RESULTS} | grep ' singletons (' | cut -d' ' -f48`)) ;
	BWA_MAPPED_OVERALL=$((${BWA_MAPPED_PP})) ;
	BWA_ALIGNED_R1=$((`echo ${BWA_RESULTS} | grep ' read1' | cut -d' ' -f 26`)) ;


	#			FILTER_ALMQUAL_RESULTS=`samtools flagstat ${TMPDIR}/bam/${TAG}.sorted.md.rel.bam ` ;
	#			FILTER_ALMQUAL_PP=$((`echo ${FILTER_ALMQUAL_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
	#			FILTER_ALMQUAL_ST=$((`echo ${FILTER_ALMQUAL_RESULTS} | grep ' singletons ' | cut -d' ' -f48`)) ;
	#			FILTER_ALMQUAL_OVERALL=$((${FILTER_ALMQUAL_ST}+${FILTER_ALMQUAL_PP})) ;
	#			FILTER_ALMQUAL_ALIGNED_R1=$((`echo ${FILTER_ALMQUAL_RESULTS} | grep ' read1' | cut -d' ' -f 26`)) ;

	ISS_RESULTS=`samtools flagstat ${TMPDIR}/bam/${TAG}.sorted.md.rel.pg.iss.bam `
	ISS_MAPPED=$((`echo ${ISS_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
	ISS_MAPPED_PP=$((`echo ${ISS_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
	ISS_MAPPED_ST=$((`echo ${ISS_RESULTS} | grep ' singletons (' | cut -d' ' -f48`)) ;
	ISS_MAPPED_OVERALL=$((${ISS_MAPPED_PP}+${ISS_MAPPED_ST})) ;
	ISS_ALIGNED_R1=$((`echo ${ISS_RESULTS} | grep ' read1' | cut -d' ' -f 26`)) ;

	ISS_FINAL=`wc -l ${TMPDIR}/bed/${TAG}.sorted.md.rel.pg.bed | cut -d' ' -f1 ` ;

	LTR_ID="${b}"
	LC_ID="${k}"

	############## DB STATS - insert - start #########################################
	echo "
	::VARIABLE SUMMARY:: dlimiters<>
	RUN_ID=<${RUN_ID}>
	RUN_NAME=<${RUN_NAME}>
	DISEASE=<${DISEASE}>
	PATIENT=<${PATIENT}>
	POOL=<${POOL}>
	TAG=<${TAG}>
	LTR_ID=<${LTR_ID}>
	LC_ID=<${LC_ID}>
	RAW_READS=<${RAW_READS}>
	QUALITY_PASSED=<${QUALITY_PASSED}>
	PHIX_MAPPING=<${PHIX_MAPPING}>
	PLASMID_MAPPED_BYPOOL=<${PLASMID_MAPPED_BYPOOL}>
	RAW_NO_PLASMID=<${RAW_NO_PLASMID}>
	BARCODE_MUX=<${BARCODE_MUX}>
	LTR_IDENTIFIED=<${LTR_IDENTIFIED}>
	TRIMMING_FINAL_LTRLC=<${TRIMMING_FINAL_LTRLC}>
	LV_MAPPED=<${LV_MAPPED}>
	BWA_MAPPED_OVERALL=<${BWA_MAPPED_OVERALL}>
	ISS_MAPPED_OVERALL=<${ISS_MAPPED_OVERALL}>
	ISS_MAPPED_PP=<${ISS_MAPPED_PP}>

	" > ${TMPDIR}/stats/${TAG}.complete.stats.txt

	echo "[TIGET] Import STATS into SUMMARY table"

done






for a in ${ASSOBCLIST[@]}; do
	TAG=${a}
	echo $TAG
	# checking VARS
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Checking Variables ";
	LV_MAPPING_RESULTS=`samtools flagstat ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${TAG}.noLTRLC.sorted.md.bam ` ; # ${TMPDIR}/bcmuxall/${DISEASE}_${PATIENT}_${POOL}.${TAG}.noLTRLC.sorted.md.bam
	LV_MAPPING_PP=$((`echo ${LV_MAPPING_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
	LV_MAPPING_ST=$((`echo ${LV_MAPPING_RESULTS} | grep ' singletons ' | cut -d' ' -f48`)) ;
	LV_MAPPING_OVERALL=$((${LV_MAPPING_ST}+${LV_MAPPING_PP})) ;
	LV_MAPPING_MAPPED=$((`echo ${LV_MAPPING_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
	LV_MAPPED=`zcat ${OUTDIR_POOL_BCMUXALL}/${DISEASE}_${PATIENT}_${POOL}.${TAG}.sorted.md.vector.list.gz | wc -l | cut -d' ' -f1 `

	PLASMID_MAPPED_BYPOOL=`zcat ${OUTDIR_VECTOR_POOL}/${DISEASE}_${PATIENT}_${POOL}.plasmids.list.gz | wc -l | cut -d' ' -f1 `

	BARCODE_MUX=$((`zcat ${TMPDIR}/bcmuxall/r1.no12.${TAG}.fastq.gz | wc -l `/4)) ;


	BWA_RESULTS=`samtools flagstat ${TMPDIR}/bam/${TAG}.sorted.md.cleaned.bam `
	BWA_INPUT_PAIRED=$((`zcat ${TMPDIR}/bcmuxall/r1.no12.${TAG}.Vector.fastq.gz  | wc -l | cut -d' ' -f1 `/4)) ;

	BWA_INPUT=$((${BWA_INPUT_PAIRED})) ;
	BWA_MAPPED=$((`echo ${BWA_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
	BWA_MAPPED_PP=$((`echo ${BWA_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
	BWA_MAPPED_ST=$((`echo ${BWA_RESULTS} | grep ' singletons (' | cut -d' ' -f48`)) ;
	BWA_MAPPED_OVERALL=$((${BWA_MAPPED_PP})) ;
	BWA_ALIGNED_R1=$((`echo ${BWA_RESULTS} | grep ' read1' | cut -d' ' -f 26`)) ;


	ISS_RESULTS=`samtools flagstat ${TMPDIR}/bam/${TAG}.sorted.md.rel.pg.iss.bam `
	ISS_MAPPED=$((`echo ${ISS_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
	ISS_MAPPED_PP=$((`echo ${ISS_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
	ISS_MAPPED_ST=$((`echo ${ISS_RESULTS} | grep ' singletons (' | cut -d' ' -f48`)) ;
	ISS_MAPPED_OVERALL=$((${ISS_MAPPED_PP}+${ISS_MAPPED_ST})) ;
	ISS_ALIGNED_R1=$((`echo ${ISS_RESULTS} | grep ' read1' | cut -d' ' -f 26`)) ;

	ISS_FINAL=`wc -l ${TMPDIR}/bed/${TAG}.sorted.md.rel.pg.bed | cut -d' ' -f1 ` ;

	STATS_TSV="${TMPDIR}/stats/summary_stats.tsv"

	# Export all needed vars to environment so we can use them in python
	export RUN_ID RUN_NAME DISEASE PATIENT POOL TAG LTR_ID LC_ID
	export RAW_READS QUALITY_PASSED PHIX_MAPPING PLASMID_MAPPED_BYPOOL
	export RAW_NO_PLASMID BARCODE_MUX LTR_IDENTIFIED TRIMMING_FINAL_LTRLC
	export LV_MAPPED BWA_MAPPED_OVERALL ISS_MAPPED_OVERALL ISS_MAPPED_PP

	python3 ${RAAVIOLIDIR}/scripts/append_stats.py "$STATS_TSV"

done
