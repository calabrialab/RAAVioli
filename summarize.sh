#!/bin/bash -x

VARIABLES_VIRAL=$1
VARIABLES_MIXED=$2
VARIABLES_REFERENCE=$3
OUTPUT_DIR=$4
BIS=$5
source $VARIABLES_VIRAL

PAR_FSAMTOOLS="772" 
PAR_FSAMTOOLS_BIS="2820"
file_par_name="${SPEC}.k${bwa_mem_k}r${bwa_mem_r}a${bwa_mem_A}t${bwa_mem_T}d${bwa_mem_d}b${bwa_mem_B}"

# ========================== aav sequences ==========================
PREFIX=".${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS}.as${bam_filter_AS}.sorted.slice.seq"
# combina tutti i file di annotazione
FIRSTROWHEADER="sample\tsourcefile\tname\tntsequence\tntsequence_len"
# init file annotation on AAV
echo -e $FIRSTROWHEADER > ${OUTPUT_DIR}/_Summary${PREFIX}.tsv

for k in $( ls ${OUTPUT_DIR}/*${PREFIX}.csv ); do
    echo $k; 
    BN=`basename $k | sed "s/${PREFIX}.csv//g"`; 
    cat $k | awk -v SAMPLE="${BN}" -v FILESRC="${k}" '{print SAMPLE"\t"FILESRC"\t"$0}' >> ${OUTPUT_DIR}/_Summary${PREFIX}.tsv
done

# ========================== aav with rearrangements ==========================
PREFIX=".${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS}.as${bam_filter_AS}.sorted.annotated.strandness"
# combina tutti i file di annotazione
FIRSTROWHEADER="sample\tsourcefile\tchr\tstart\tend\tname\tscore\tstrand\tcigar\taav_element_chr\taav_element_start\taav_element_end\taav_element_name\taav_element_score\taav_element_strand\taav_element_idname\taav_element_distance"
# init file annotation on AAV
echo -e $FIRSTROWHEADER > ${OUTPUT_DIR}/_Summary${PREFIX}.tsv


for k in $( ls ${OUTPUT_DIR}/*${PREFIX}.bed ); do
    echo $k; 
    BN=`basename $k | sed "s/${PREFIX}.bed//g"`; 
    cat $k | awk -v SAMPLE="${BN}" -v FILESRC="${k}" '{print SAMPLE"\t"FILESRC"\t"$0}' >> ${OUTPUT_DIR}/_Summary${PREFIX}.tsv
done

# ========================== aav with rearrangements and all features annotated ==========================
PREFIX=".${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS}.as${bam_filter_AS}.sorted.annotated.allfeatures.strandness"
# combina tutti i file di annotazione
FIRSTROWHEADER="sample\tsourcefile\tchr\tstart\tend\tname\tscore\tstrand\tcigar\taav_element_chr\taav_element_start\taav_element_end\taav_element_name\taav_element_score\taav_element_strand\taav_element_idname\taav_element_distance"
# init file annotation on AAV
echo -e $FIRSTROWHEADER > ${OUTPUT_DIR}/_Summary${PREFIX}.tsv


for k in $( ls ${OUTPUT_DIR}/*${PREFIX}.bed ); do
echo $k; 
BN=`basename $k | sed "s/${PREFIX}.bed//g"`; 
cat $k | awk -v SAMPLE="${BN}" -v FILESRC="${k}" '{print SAMPLE"\t"FILESRC"\t"$0}' >> ${OUTPUT_DIR}/_Summary${PREFIX}.tsv
done

# ========================== aav unique best ==========================
PREFIX=".${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS_BIS}.as${bam_filter_AS}.sorted.annotated.strandness"
# combina tutti i file di annotazione
FIRSTROWHEADER="sample\tsourcefile\tchr\tstart\tend\tname\tscore\tstrand\tcigar\taav_element_chr\taav_element_start\taav_element_end\taav_element_name\taav_element_score\taav_element_strand\taav_element_idname\taav_element_distance"
# init file annotation on AAV
echo -e $FIRSTROWHEADER > ${OUTPUT_DIR}/_Summary${PREFIX}.tsv


for k in $( ls ${OUTPUT_DIR}/*${PREFIX}.bed ); do
    echo $k; 
    BN=`basename $k | sed "s/${PREFIX}.bed//g"`; 
    cat $k | awk -v SAMPLE="${BN}" -v FILESRC="${k}" '{print SAMPLE"\t"FILESRC"\t"$0}' >> ${OUTPUT_DIR}/_Summary${PREFIX}.tsv
done

# ========================== aav unique best coverage ==========================
PREFIX=".${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS}.as${bam_filter_AS}.sorted.coverageByElement"
# combina tutti i file di annotazione
FIRSTROWHEADER="sample\tsourcefile\taav_element_chr\taav_element_start\taav_element_end\taav_element_name\taav_element_score\taav_element_strand\taav_element_idname\tn_covered_reads\tn_covered_bases\taav_element_size\tcovered_bases_on_element_size"
# init file annotation on AAV
echo -e $FIRSTROWHEADER > ${OUTPUT_DIR}/_Summary${PREFIX}.tsv


for k in $( ls ${OUTPUT_DIR}/*${PREFIX}.bed ); do
echo $k; 
BN=`basename $k | sed "s/${PREFIX}.bed//g"`; 
NREADS=`basename $k | sed "s/${PREFIX}.bed//g"`; 
cat $k | awk -v SAMPLE="${BN}" -v FILESRC="${k}" '{print SAMPLE"\t"FILESRC"\t"$0}' >> ${OUTPUT_DIR}/_Summary${PREFIX}.tsv
done



file_par_name="${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS}.as${bam_filter_AS}"
if  [ "$#" -eq 5 ]
then
    # ========================== reference ==========================
    # now on reference genomes

    source $VARIABLES_REFERENCE

    second_file_par="slice.${SPEC}.q${sam_view_q}.${PAR_FSAMTOOLS_BIS}.as${bam_filter_AS}.nm${bam_filter_NM}"


    PREFIX=".${file_par_name}.sorted.${second_file_par}.sorted.annotated"
    FIRSTROWHEADER="sample\tsourcefile\tchr\tstart\tend\tname\tscore\tstrand\tcigar\tgene_chr\tgene_annotationsource\tgene_elementtype\tgene_start\tgene_end\tgene_f1\tgene_strand\tgene_f2\tgene_details\tdistance_to_gene"
    echo -e $FIRSTROWHEADER > ${OUTPUT_DIR}/_Summary${PREFIX}.tsv

    for k in $( ls ${OUTPUT_DIR}/*${PREFIX}.bed ); do
    echo $k; 
    BN=`basename $k | sed "s/${PREFIX}.bed//g"`; 
    cat $k | awk -v SAMPLE="${BN}" -v FILESRC="${k}" '{print SAMPLE"\t"FILESRC"\t"$0}' >> ${OUTPUT_DIR}/_Summary${PREFIX}.tsv
    done



    # ========================== reference ========================== full gene ==========================
    # now on reference genomes
    PREFIX=".${file_par_name}.sorted.${second_file_par}.sorted.annotatedFullGene.strandness"
    FIRSTROWHEADER="sample\tsourcefile\tchr\tstart\tend\tname\tscore\tstrand\tcigar\tgene_chr\tgene_start\tgene_end\tgene_name2\tgene_score\tgene_strand\tdistance_to_gene"
    echo -e $FIRSTROWHEADER > ${OUTPUT_DIR}/_Summary${PREFIX}.tsv

    for k in $( ls ${OUTPUT_DIR}/*${PREFIX}.bed ); do
    echo $k; 
        BN=`basename $k | sed "s/${PREFIX}.bed//g"`; 
        cat $k | awk -v SAMPLE="${BN}" -v FILESRC="${k}" '{print SAMPLE"\t"FILESRC"\t"$0}' >> ${OUTPUT_DIR}/_Summary${PREFIX}.tsv
    done

    # ========================== reference ========================== strandness ========================== 
    # now on reference genomes
    PREFIX=".${file_par_name}.sorted.${second_file_par}.sorted.annotated.strandness"
    FIRSTROWHEADER="sample\tsourcefile\tchr\tstart\tend\tname\tscore\tstrand\tcigar\tgene_chr\tgene_annotationsource\tgene_elementtype\tgene_start\tgene_end\tgene_f1\tgene_strand\tgene_f2\tgene_details\tdistance_to_gene"
    echo -e $FIRSTROWHEADER > ${OUTPUT_DIR}/_Summary${PREFIX}.tsv

    for k in $( ls ${OUTPUT_DIR}/*${PREFIX}.bed ); do
    echo $k; 
    BN=`basename $k | sed "s/${PREFIX}.bed//g"`; 
    cat $k | awk -v SAMPLE="${BN}" -v FILESRC="${k}" '{print SAMPLE"\t"FILESRC"\t"$0}' >> ${OUTPUT_DIR}/_Summary${PREFIX}.tsv
    done
fi

# ========================== reference mixed aav ========================== strandness ========================== 
# now on reference genomes

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


