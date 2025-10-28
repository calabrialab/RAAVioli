counter=0
MAXTHREADS=1;
if [ -z "$CHRV_NAME" ]
then
  CHRV_NAME="chrV"
fi

if [ -z "$MAXCLUSTERD" ]
then
  MAXCLUSTERD=20
fi

if [ -z "$MINAAVMATCHES" ]
then
  MINAAVMATCHES=30
fi

if [ -z "$MERGECOL" ]
then
  MERGECOL="CompleteAmplificationID"
fi

#if [ -z "$POOLSAID" ]
#then
#  POOLSAID=${POOL}
#fi

for TAG in ${ASSOBCLIST[@]}; do
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET]  Searching for chimeras"
  python3 ${RAAVIOLIDIR}/scripts/adaptive.v11.py -i ${TMPDIR}/bam/${POOLSAID}/${TAG}.sorted.md.rel.pg.iss.bam -F ${TMPDIR}/bam/${POOLSAID}/${TAG}.F4.sorted.md.bam -o ${TMPDIR}/iss/${POOLSAID}/${TAG}.results -l 50 -g 50 -s ${SUBOPTH} -c ${CHRV_NAME} &
	counter=$((counter + 1))
	if (( counter >= MAXTHREADS )); then
		wait
		counter=0
	fi
done
wait

if [ -z "$OUTPUT_NAME" ]
then
  OUTPUT_NAME="v1"
fi
output_fn="results.${OUTPUT_NAME}.joined.tsv.gz";
python3 ${RAAVIOLIDIR}/scripts/results_concat.py -i ${TMPDIR}/iss/${POOLSAID}/ -o ${TMPDIR}/iss/${POOLSAID}/${output_fn} -a ${ASSOCIATIONFILE}
rm  ${TMPDIR}/iss/${POOLSAID}/*.results.R1.tsv

output_fn_mod=$(echo "$output_fn" | sed 's/.tsv.gz/.modif.tsv.gz/')
final_output_name=$output_fn_mod
echo $output_fn_mod
if [ -z "$ITR_DF" ]
then
    echo "modifying"
    python3 ${RAAVIOLIDIR}/scripts/0.1.modify_raw_matrix.py -i ${TMPDIR}/iss/${POOLSAID}/${output_fn} -o ${TMPDIR}/iss/${POOLSAID}/${output_fn_mod} -a ${ASSOCIATIONFILE}  -n ${POOL} -m ${MINAAVMATCHES}
else
    echo "modfifying"
    python3 ${RAAVIOLIDIR}/scripts/0.1.modify_raw_matrix.py -i ${TMPDIR}/iss/${POOLSAID}/${output_fn} -o ${TMPDIR}/iss/${POOLSAID}/${output_fn_mod} -a ${ASSOCIATIONFILE}  -n ${POOL} -m ${MINAAVMATCHES} -d ${ITR_DF}
fi

if  [ -n "${system_sequences_df}" ]
then
  mkdir ${TMPDIR}/iss/${POOLSAID}/expo_analysis
  python3 ${RAAVIOLIDIR}/scripts/0.2.create_fasta_dividedbyadapters.py -i ${TMPDIR}/iss/${POOLSAID}/${output_fn_mod} -o ${TMPDIR}/iss/${POOLSAID}/expo_analysis/  -s ${system_sequences_df}
  for file in ${TMPDIR}/iss/${POOLSAID}/expo_analysis/*_expo_reads.fa;
  do
      bn=`basename ${file} | sed 's/_expo_reads.fa//g'`
      echo ${file} ${bn}
     flexbar -r ${file} -b ${TMPDIR}/iss/${POOLSAID}/expo_analysis/${bn}.fa --barcode-trim-end LEFT  -be 0.05 -m 2 -O ${TMPDIR}/iss/${POOLSAID}/expo_analysis/${bn}.log -t ${TMPDIR}/iss/${POOLSAID}/expo_analysis/flexbarOut

  done
  output_fn_mod_check=$(echo "$output_fn_mod" | sed 's/.tsv.gz/.checked.tsv.gz/')
  echo $output_fn_mod_check
  python3 ${RAAVIOLIDIR}/scripts/0.3.remove_expois.py -i ${TMPDIR}/iss/${POOLSAID}/${output_fn_mod} -o ${TMPDIR}/iss/${POOLSAID}/${output_fn_mod_check} -f ${TMPDIR}/iss/${POOLSAID}/expo_analysis/
  final_output_name=$output_fn_mod_check
fi

mkdir ${TMPDIR}/matrix
mkdir ${TMPDIR}/matrix/${POOLSAID}
python3 ${RAAVIOLIDIR}/scripts/1.create_matrix_clustering_sorted.py -i ${TMPDIR}/iss/${POOLSAID}/${final_output_name} -o ${TMPDIR}/matrix/${POOLSAID}/ -b ${OUTPUT_NAME} -c ${MAXCLUSTERD} -m ${MERGECOL}
#python3 ${RAAVIOLIDIR}/scripts/1.create_matrix_clustering_sorted_onlyIS.py -i ${TMPDIR}/iss/${POOLSAID}/${final_output_name} -o ${TMPDIR}/matrix/${POOLSAID}/ -b ${OUTPUT_NAME}
if [ -n "${ANNOTATIONGTF}" ]
then
  bash ${RAAVIOLIDIR}/utils_scripts/annotate_matrix.sh  -m ${TMPDIR}/matrix/${POOLSAID}/SeqCount_${OUTPUT_NAME}.CLUSTER${MAXCLUSTERD}.tsv -g ${ANNOTATIONGTF} -t vispa2 -o ${TMPDIR}/matrix/${POOLSAID}/
  bash ${RAAVIOLIDIR}/utils_scripts/annotate_matrix.sh  -m ${TMPDIR}/matrix/${POOLSAID}/ShsCount_${OUTPUT_NAME}.CLUSTER${MAXCLUSTERD}.tsv -g ${ANNOTATIONGTF} -t vispa2 -o ${TMPDIR}/matrix/${POOLSAID}/
  bash ${RAAVIOLIDIR}/utils_scripts/annotate_matrix.sh  -m ${TMPDIR}/matrix/${POOLSAID}/SeqCount_${OUTPUT_NAME}.CLUSTER${MAXCLUSTERD}.cleaned1.tsv -g ${ANNOTATIONGTF} -t vispa2 -o ${TMPDIR}/matrix/${POOLSAID}/
  bash ${RAAVIOLIDIR}/utils_scripts/annotate_matrix.sh  -m ${TMPDIR}/matrix/${POOLSAID}/ShsCount_${OUTPUT_NAME}.CLUSTER${MAXCLUSTERD}.cleaned1.tsv -g ${ANNOTATIONGTF} -t vispa2 -o ${TMPDIR}/matrix/${POOLSAID}/
  cp ${TMPDIR}/matrix/${POOLSAID}/*.annotated.tsv.gz ${OUTDIR_POOL_MATRIX}
fi
cp ${TMPDIR}/iss/${POOLSAID}/${final_output_name}  ${OUTDIR_POOL_MATRIX}
cp ${TMPDIR}/matrix/${POOLSAID}/*.tsv ${OUTDIR_POOL_MATRIX}

