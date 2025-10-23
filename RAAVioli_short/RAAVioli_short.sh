RAAVIOLIDIR=$(pwd)
echo $RAAVIOLIDIR
export RAAVIOLIDIR
# Activate conda environment
ENV_NAME="RAAVioliShort_env"
if command -v conda >/dev/null 2>&1; then
    # Load conda shell integration
    eval "$(conda shell.bash hook)"
    conda activate "${ENV_NAME}" || {
        echo "[ERROR] Failed to activate conda environment: ${ENV_NAME}"
        exit 1
    }
else
    echo "[ERROR] conda not found. Please install Miniconda/Mamba and run setup first."
    exit 1
fi

# Check if the file is provided as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <mandatoryvars.txt>"
    exit 1
fi

mandatoryvars=$1

source $mandatoryvars


OUTDIR_MERGE_MATRIX="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/matrix/";
OUTDIR_POOL_MATRIX="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/matrix/"${POOL};
OUTDIR_MERGE_STATS="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/stats/";
OUTDIR_POOL_STATS="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/stats/"${POOL};

mkdir ${NGSWORKINGPATH}
mkdir ${NGSWORKINGPATH}/${DISEASE}
mkdir ${NGSWORKINGPATH}/${DISEASE}/${PATIENT}
mkdir ${OUTDIR_MERGE_MATRIX};
mkdir ${OUTDIR_POOL_MATRIX};
mkdir ${OUTDIR_MERGE_STATS};
mkdir ${OUTDIR_POOL_STATS};


mkdir ${TMPDIR} ;
TMPDIR="${TMPDIR}/${POOL}"
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

#BCLTRf=`cut "${BARCODE_LTR}" -f1`;
#BCLTR=(${BCLTRf}) ;
#BCLCf=`cut "${BARCODE_LC}" -f1`;
#BCLC=(${BCLCf}) ;



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

ASSOBCLIST=(`cut -f${target_index} ${ASSOCIATIONFILE} | tail -n+2 | sort | uniq `);



##### ================ COPY DATA INTO TMP DIR ================== #####

if [ -n "$alignment_vars_file" ]; then
	source $alignment_vars_file

		for TAG in ${ASSOBCLIST[@]}; do
		##### ==========================  MOVING READS IN THE BCMUXALL DIR IF STEP 1 AND 2 NOT DONE ================================= #####
			echo "MOVING READS IN THE BCMUXALL DIR"
			cp ${R1_FASTQ}/${TAG}.r1.fastq.gz  ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_1.fastq.gz
			cp ${R2_FASTQ}/${TAG}.r2.fastq.gz  ${TMPDIR}/bcmuxall/${TAG}.cleanR1LCR2LC_2.fastq.gz
    	done
	source step_alignment.sh
fi



if [ -n "$isr_vars_file" ]; then
  source $isr_vars_file
  source step_is_identification.sh
fi

#### clean tmp dir
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Cleaning TMP Directory... "
if [ ${REMOVE_TMP_DIR} = "remove_tmp_yes" ]
	then
	rm -fr ${TMPDIR};
fi

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Completed ";

