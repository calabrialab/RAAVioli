#!/bin/bash -x

VARIABLES_VIRAL=$1
VARIABLES_MIXED=$2
OUTPUT_DIR=$3
VARIABLES_STEPR=$4
INPUT_FILE_TSV=$5
source $VARIABLES_VIRAL

PAR_FSAMTOOLS="772" 
file_par_name="${SPEC}.k${bwa_mem_k}r${bwa_mem_r}a${bwa_mem_A}t${bwa_mem_T}d${bwa_mem_d}b${bwa_mem_B}"
file_par_name="${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS}.as${bam_filter_AS}"



# ========================== reference mixed aav ========================== strandness ========================== 

source $VARIABLES_MIXED
second_file_par="slice.${SPEC}.q${sam_view_q}.${PAR_FSAMTOOLS}.as${bam_filter_AS}.nm${bam_filter_NM}"

PREFIX=".${file_par_name}.sorted.${second_file_par}.sorted.annotated.strandness"
FIRSTROWHEADER="sample\tsourcefile\tchr\tstart\tend\tname\tscore\tstrand\tcigar\tgene_chr\tgene_annotationsource\tgene_elementtype\tgene_start\tgene_end\tgene_f1\tgene_strand\tgene_f2\tgene_details\tdistance_to_gene"
echo -e $FIRSTROWHEADER > ${OUTPUT_DIR}/_Summary${PREFIX}.tsv

for k in $( ls ${OUTPUT_DIR}/*${PREFIX}.bed ); do
    echo $k; 
    BN=`basename $k | sed "s/${PREFIX}.bed//g"`; 
    cat $k | awk -v SAMPLE="${BN}" -v FILESRC="${k}" '{print SAMPLE"\t"FILESRC"\t"$0}' >> ${OUTPUT_DIR}/_Summary${PREFIX}.tsv
done

source $VARIABLES_STEPR
mkdir -p ${OUTPUT_DIR}/final_results
PIPEDIR=$(pwd)
Rscript step_R.R -o ${OUTPUT_DIR} -c ${PIPEDIR} -i ${INPUT_FILE_TSV} -s ${OUTPUT_DIR}/_Summary${PREFIX}.tsv -m $min_cigar_alm_width -n $reads_all_alm_by_cigar_chimera_details_alm_size

