#!/bin/bash
source config.txt
source $6

INPUT_FILE=$1
file_par_name=$2
MAXTHREADS=$3
PAR_FSAMTOOLS=$4
GENOME=$5
OUTPUT_DIR=$7
ANNOTATION=$8

fq_files=($(awk -F$'\t' 'NR>=2 {print $NF}' ${INPUT_FILE}))
list_bn=()

echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Align to ${SPEC} genome ============"
for fq_file in "${fq_files[@]}"
do
    BN=`basename $fq_file | sed 's/.fastq.gz//g'`; 
    BN=${BN}.${file_par_name}.sorted
    list_bn+=($BN)
    bwa_mem_R="@RG\tID:${BN}\tSM:${SPEC}.k${bwa_mem_k}r${bwa_mem_r}a${bwa_mem_A}t${PAR_T}b${PAR_B}\tCN:TIGET"
    
    $BWA  mem -k ${bwa_mem_k} -r ${bwa_mem_r} -A ${bwa_mem_A} -T ${bwa_mem_T} -d ${bwa_mem_d} -B ${bwa_mem_B} -O ${bwa_mem_O} \
     -E ${bwa_mem_E} -L ${bwa_mem_L} -R ${bwa_mem_R} -x ${bwa_mem_x} -t ${MAXTHREADS} ${GENOME} <( zcat ${OUTPUT_DIR}/${BN}.slice.fastq.gz ) | \
      $SAMTOOLS view -F ${PAR_FSAMTOOLS} -q ${sam_view_q} -uS - | \
       $SAMTOOLS sort - -o ${OUTPUT_DIR}/${BN}.slice.${SPEC}.q${sam_view_q}.${PAR_FSAMTOOLS}.sorted.bam 
done


second_file_par="slice.${SPEC}.q${sam_view_q}.${PAR_FSAMTOOLS}"

echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering ============"

for BN in "${list_bn[@]}"
do 
    $BAMTOOLS  filter -tag "AS:>=${bam_filter_AS}" -tag "NM:<=${bam_filter_NM}" -in ${OUTPUT_DIR}/${BN}.${second_file_par}.sorted.bam \
     -out ${OUTPUT_DIR}/${BN}.${second_file_par}.as${bam_filter_AS}.nm${bam_filter_NM}.sorted.bam
    
    $SAMTOOLS index ${OUTPUT_DIR}/${BN}.${second_file_par}.as${bam_filter_AS}.nm${bam_filter_NM}.sorted.bam &
done
wait

#echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Create BigWig file ============"
#for BN in "${list_bn[@]}"
#do
#    $BAMCOVERAGE    -b ${OUTPUT_DIR}/${BN}.${second_file_par}.as${bam_filter_AS}.nm${bam_filter_NM}.sorted.bam -o \
#     ${OUTPUT_DIR}/${BN}.${second_file_par}.as${bam_filter_AS}.nm${bam_filter_NM}.sorted.bw &
#done
#wait

echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Create BED file ============"
for BN in "${list_bn[@]}"
do
    $BEDTOOLS   bamtobed -cigar -i ${OUTPUT_DIR}/${BN}.${second_file_par}.as${bam_filter_AS}.nm${bam_filter_NM}.sorted.bam > \
     ${OUTPUT_DIR}/${BN}.${second_file_par}.as${bam_filter_AS}.nm${bam_filter_NM}.bed 
     sort -k1,1 -k2,2n -k3,3n ${OUTPUT_DIR}/${BN}.${second_file_par}.as${bam_filter_AS}.nm${bam_filter_NM}.bed > ${OUTPUT_DIR}/${BN}.${second_file_par}.as${bam_filter_AS}.nm${bam_filter_NM}.sorted.bed 
done


echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Annotate BED file with closest gene (our annotation) ============"
for BN in "${list_bn[@]}"
do

    $BEDTOOLS   closest -b ${ANNOTATION} -a ${OUTPUT_DIR}/${BN}.${second_file_par}.as${bam_filter_AS}.nm${bam_filter_NM}.sorted.bed \
     -d -D ref -t first > ${OUTPUT_DIR}/${BN}.${second_file_par}.as${bam_filter_AS}.nm${bam_filter_NM}.sorted.annotated.strandness.bed 

done

