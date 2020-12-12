#!/bin/bash

source config.txt
source $VARIABLES_VIRAL
helpFunction()
{
   echo ""
   echo "Sample usage: $0 -i sample_label.tsv -t threads -v viral_genome.fa -r reference.fa -R 1 -A annotation.bed -o output_dir -m mixed_genome.fa"
   echo -e "\n"
   echo -e "\t-i the .tsv file with the paths to fastq files as last column.\n"
   echo -e "\t-t max threads to be used.\n"
   echo -e "\t-v (optional) the fasta file with viral genome (e.g. AAV).\n\t   The bwa-index will be created in the same directory. \n\t   If you have already an index please see -V.\n\t   You must specify -V if you don't specify -v.\n"
   echo -e "\t-r (optional) the fasta file with the reference genome (e.g. hg19).\n\t   The bwa-index will be created in the same directory. \n\t   If you have already a bwa-index please see -R. N.B.\n\t   You must specify -R if you don't specify -r.\n"
   echo -e "\t-V (optional) path to the viral bwa-index with basename\n\t   (e.g. if you have the index in /home/resources/genome/index\n\t   directory and it has as basename aav.fa\n\t   you have to specify home/resources/genome/index/aav.fa ).\n\t   If specified the index of the viral genome will not be made.\n\t   If you don't specify -V you must specify -v.\n"
   echo -e "\t-R (optional) path to the reference bwa-index with basename\n\t   (e.g. if you have the index in /home/resources/genome/index\n\t   directory and it has as basename hg19.fa\n\t   you have to specify home/resources/genome/index/hg19.fa ).\n\t   If specified the index of the reference genome will not be made.\n\t   If you don't specify -R you must specify -r.\n"
   echo -e "\t-m (optional) the fasta file with the mixed genome\n\t   N.B. viral genome must be appended at the end of reference genome\n\t   with the sequence name chrV.\n\t   Please note that if not specified it will be created and \n\t   you must specify -v and -r \n\t   (since index could be located in a different dir\n\t   and to create the mixed genome both genomes are needed). \n\t   In this case if you already have\n\t   the viral index and/or the reference index \n\t   in the same directory you can specify -V 1 and/or -R 1 instead \n\t   of specifying twice the same path for -v and -V (or -r and -R).\n"
   echo -e "\t-M (optional) bwa-index of the mixed_genome.\n\t   If specified you can omit -m.\n" 
   echo -e "\t-a the bed file with the custom annotation.\n"
   echo -e "\t-o path to the output directory.\n"
   echo -e "\t-b (optional) any value. If specified also 2820 will be made.\n"
   echo -e "\t Please read the Read.me to have more detailed info.\n"
   exit 1 
}
while getopts "i:t:v:r:m:a:o:V:R:M:b:" opt
do
   case "$opt" in
      i ) INPUT_FILE="$OPTARG" ;;
      t ) MAXTHREADS="$OPTARG" ;;
      v ) VIRALGENOME="$OPTARG" ;;
      r ) REFGENOME="$OPTARG" ;;
      a ) ANNOTATIONGTF="$OPTARG" ;;
      o ) OUTPUT_DIR="$OPTARG" ;;
      m ) MIXEDGENOME="$OPTARG" ;;
      V ) VIRALINDEX="$OPTARG" ;;
      R ) REFINDEX="$OPTARG" ;;
      M ) MIXEDINDEX="$OPTARG" ;;
      b ) BIS="$OPTARG" ;;
      ? ) helpFunction ;;
   esac
done


if [ -z "$INPUT_FILE" ] || [ -z "$MAXTHREADS" ] || [ -z "$ANNOTATIONGTF" ] || [ -z "$OUTPUT_DIR" ] 
then
   echo "Parameters missing!";
   helpFunction
fi


if [ -z "$VIRALGENOME" ] && [ -z "$VIRALINDEX" ] 
then
   echo "You must specify -v or -V or both!"
   helpFunction
fi

if [ -z "$REFGENOME" ] && [ -z "$REFINDEX" ] 
then
   echo "You must specify -r or -R or both!"
   helpFunction
fi

# checking if we are able to create the mixed genome if not specified
if [[ (-z "$MIXEDGENOME" && -z "$MIXEDINDEX") && (-z "$VIRALGENOME" || -z "$REFGENOME") ]] 
then
    echo "You must specify -v and -r if neither -m nor -M are specified."
    exit 1
fi


# checking if files exist
if [ ! -s "$INPUT_FILE" ]
then
   echo "${INPUT_FILE} does not exist or has size zero"
   exit 1
fi
if [ "$VIRALINDEX" = "1" ] && [ -z "$VIRALGENOME" ]
then
    echo "You specified -V 1 but not -v. Specifying -V 1 means that the index is in the same location of the genome specified in -v."
    helpFunction
fi
if [ "$REFINDEX" = "1" ] && [ -z "$REFGENOME" ]
then
    echo "You specified -R 1 but not -r. Specifying -R 1 means that the index is in the same location of the genome specified in -r."
    helpFunction
fi

if [ ! -z "$VIRALGENOME" ] && [ ! -s "$VIRALGENOME" ]
then
   echo "${VIRALGENOME} does not exist or has size zero"
   exit 1
fi

if [ ! -z "$VIRALGENOME" ]
then
    var=$(grep ">" ${VIRALGENOME})
    var=`echo $var | sed 's/ *$//g'`
    if [ "$var" != ">chrV" ]
    then
        echo "${VIRALGENOME} must be a single sequence with sequence name equal to >chrV"
        exit 1
    fi 
fi
if [ ! -z "$REFGENOME" ] && [ ! -s "$REFGENOME" ]
then
   echo "${REFGENOME} does not exist or has size zero"
   exit 1
fi


if [ ! -d "$OUTPUT_DIR" ]
then
    mkdir $OUTPUT_DIR
    mkdir $OUTPUT_DIR/resources
elif [ ! -d "$OUTPUT_DIR/resources" ]
then
    mkdir $OUTPUT_DIR/resources
fi

# If mixed index is specified we don't need mixed genome. If neither MIXEDINDEX nor MIXEDGENOME is specified we have to create both
if [ ! -z "$MIXEDINDEX" ]
then
    MIXEDGENOME=$MIXEDINDEX
else
    if [ ! -z "$MIXEDGENOME" ]
    then
        $BWA index -a bwtsw ${MIXEDGENOME}
    else
        echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating mixed genome in ${OUTPUT_DIR}/resources ============"
        cp $REFGENOME $OUTPUT_DIR/resources/mixed.fa
        MIXEDGENOME="$OUTPUT_DIR/resources/mixed.fa"
        cat ${VIRALGENOME} >>  ${MIXEDGENOME}
        echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Indexing mixed genome ============"
        $BWA index -a bwtsw ${MIXEDGENOME}
        echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Done ============"
    fi
fi

if [ -z "$VIRALINDEX" ]
then
    echo "Indexing ${VIRALGENOME}"
    #$BWA index -a bwtsw ${VIRALGENOME}
elif [ -z "$VIRALGENOME" ] 
then
    VIRALGENOME=$VIRALINDEX
elif [ "$VIRALINDEX" != "1" ]
then
    VIRALGENOME=$VIRALINDEX
fi

    

PAR_FSAMTOOLS="772" 
PAR_FSAMTOOLS_BIS="2820"


#reading all paths from .tsv file

fq_files=($(awk -F$'\t' 'NR>=2 {print $NF}' ${INPUT_FILE}))
file_par_name="${SPEC}.k${bwa_mem_k}r${bwa_mem_r}a${bwa_mem_A}t${bwa_mem_T}d${bwa_mem_d}b${bwa_mem_B}"
list_bn=()

echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Align to Vector genome and Filtering ============"
for fq_file in "${fq_files[@]}"
do
    BN=`basename $fq_file | sed 's/.fastq.gz//g'`; 
    bwa_mem_R="@RG\tID:${BN}${bwa_mem_R}\tCN:TIGET"
    list_bn+=($BN)
    
    $BWA mem -k ${bwa_mem_k} -r ${bwa_mem_r} -A ${bwa_mem_A} -T ${bwa_mem_T} -d ${bwa_mem_d} -B ${bwa_mem_B} -O ${bwa_mem_O} \
     -E ${bwa_mem_E} -L ${bwa_mem_L} -R ${bwa_mem_R} -t ${MAXTHREADS} ${VIRALGENOME} <( zcat ${fq_file} ) | \
       $SAMTOOLS view -F ${PAR_FSAMTOOLS} -q $sam_view_q -uS - | \
        $SAMTOOLS sort - -o ${OUTPUT_DIR}/${BN}.${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS}.sorted.bam &   
    
done
wait

file_par_name="${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS}"

for BN in "${list_bn[@]}"
do
    $BAMTOOLS  filter -tag "AS:>=${bam_filter_AS}" -in ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bam -out \
    ${OUTPUT_DIR}/${BN}.${file_par_name}.as${bam_filter_AS}.sorted.bam &
done
wait
file_par_name="${file_par_name}.as${bam_filter_AS}"
for BN in "${list_bn[@]}"
do
    $SAMTOOLS index ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bam &
done
wait

echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Create BED file ============"
for BN in "${list_bn[@]}"
do
    $BEDTOOLS   bamtobed -cigar -i ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bam > ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bed &
done
wait



echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Annotate BED file with closest gene (our annotation) ============"
for BN in "${list_bn[@]}"
do
    $BEDTOOLS   closest -b ${ANNOTATIONGTF} -a ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bed -s -d -D ref -t first > \
    ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.annotated.bed &
    $BEDTOOLS   closest -b ${ANNOTATIONGTF} -a ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bed -d -D ref -t first > \
    ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.annotated.strandness.bed &
    $BEDTOOLS   closest -b ${ANNOTATIONGTF} -a ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bed -d -D ref > \
    ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.annotated.allfeatures.strandness.bed &
    $BEDTOOLS   coverage -b ${ANNOTATIONGTF} -a ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bed > \
    ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.coverageByElement.bed &
    wait
done


#create 2820 from 772
#Doing everything again with a different PAR_FSAMTOOLS
file_par_name="${SPEC}.k${bwa_mem_k}r${bwa_mem_r}a${bwa_mem_A}t${bwa_mem_T}d${bwa_mem_d}b${bwa_mem_B}"


echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering ============"
for BN in "${list_bn[@]}"
do
    $SAMTOOLS view ${OUTPUT_DIR}/${BN}.${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS}.sorted.bam -F ${PAR_FSAMTOOLS_BIS} -q ${sam_view_q} -uS |\
    $SAMTOOLS sort - -o ${OUTPUT_DIR}/${BN}.${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS_BIS}.sorted.bam &   
done
wait
PAR_FSAMTOOLS=${PAR_FSAMTOOLS_BIS}
file_par_name="${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS}"

for BN in "${list_bn[@]}"
do
    $BAMTOOLS  filter -tag "AS:>=${bam_filter_AS}" -in ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bam -out \
    ${OUTPUT_DIR}/${BN}.${file_par_name}.as${bam_filter_AS}.sorted.bam &
done
wait
file_par_name="${file_par_name}.as${bam_filter_AS}"
for BN in "${list_bn[@]}"
do
    $SAMTOOLS index ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bam &
done
wait

echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Create BED file ============"
for BN in "${list_bn[@]}"
do
    $BEDTOOLS   bamtobed -cigar -i ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bam > ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bed &
done
wait



echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Annotate BED file with closest gene (our annotation) ============"
for BN in "${list_bn[@]}"
do
    $BEDTOOLS   closest -b ${ANNOTATIONGTF} -a ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bed -s -d -D ref -t first > \
    ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.annotated.bed &
    $BEDTOOLS   closest -b ${ANNOTATIONGTF} -a ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bed -d -D ref -t first > \
    ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.annotated.strandness.bed &
    $BEDTOOLS   closest -b ${ANNOTATIONGTF} -a ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bed -d -D ref > \
    ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.annotated.allfeatures.strandness.bed &
    $BEDTOOLS   coverage -b ${ANNOTATIONGTF} -a ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bed > \
    ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.coverageByElement.bed &
    wait
done



PAR_FSAMTOOLS="772" 
file_par_name="${SPEC}.k${bwa_mem_k}r${bwa_mem_r}a${bwa_mem_A}t${bwa_mem_T}d${bwa_mem_d}b${bwa_mem_B}"
file_par_name="${file_par_name}.q${sam_view_q}F${PAR_FSAMTOOLS}"
file_par_name="${file_par_name}.as${bam_filter_AS}"


for fq_file in "${fq_files[@]}"
do
    BN=`basename $fq_file | sed 's/.fastq.gz//g'`; 
    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Extract reads from raw data"
    cat ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.bed | cut -f4 | sort | uniq > ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.headerlist
    zcat ${fq_file} | python $FQEXTRACT ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.headerlist | \
    pigz -f -c > ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.slice.fastq.gz 

    echo "[AP] ============ <`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Get sequence file ============"
    zcat ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.slice.fastq.gz | $FASTQ_TO_FASTA -Q33 | ruby $FASTA_TO_CSV | tr " " "\t" | \
    awk '{ print $0"\t"length($2) }' > ${OUTPUT_DIR}/${BN}.${file_par_name}.sorted.slice.seq.csv
done


bash step2.sh ${INPUT_FILE} ${file_par_name} ${MAXTHREADS} ${PAR_FSAMTOOLS} ${MIXEDGENOME} ${VARIABLES_MIXED} ${OUTPUT_DIR} 

#checking if BIS is setted. If yes we check if we have to index reference genome (-rIndex 1) and then call the step2 for reference genome.
if [ ! -z "$BIS" ]
then
    if [ -z "$REFINDEX" ]
        then
            echo "Indexing ${REFGENOME}"
            $BWA index -a bwtsw ${REFGENOME}
        elif [ -z "$REFGENOME" ] 
        then
            REFGENOME=$REFINDEX
        elif [ "$REFINDEX" != "1" ]
        then
            REFGENOME=$REFINDEX
    fi
    bash step2.sh ${INPUT_FILE} ${file_par_name} ${MAXTHREADS} ${PAR_FSAMTOOLS_BIS} ${REFGENOME} ${VARIABLES_REFERENCE} ${OUTPUT_DIR} 
fi


bash summarize.sh $VARIABLES_VIRAL $VARIABLES_MIXED $VARIABLES_REFERENCE $OUTPUT_DIR $BIS
