
# FASTA of RANDOM TAGs
if [ ${FASTQ_QF} = "slim" ]; then
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Export Random TAGs from R2"
	trimmomatic SE -phred33 -threads $MAXTHREADS ${R2_FASTQ} "/dev/stdout" CROP:12 | fastq_to_fasta -Q 33 -n | pigz -f -c > ${OUTDIR_POOL_QUAL}/r2.${POOL}.qf.noPlasmids.noPhiX.TAGs.fa.gz
fi

## remove first N bases
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Remove first 12 bases"
trimmomatic SE -threads ${MAXTHREADS} -phred33 ${R1_FASTQ} ${TMPDIR}/r1.no12.fastq.gz HEADCROP:12 &
trimmomatic SE -threads ${MAXTHREADS} -phred33 ${R2_FASTQ} ${TMPDIR}/r2.no12.fastq.gz HEADCROP:12 &
wait

#rm ${R1_FASTQ} ${R2_FASTQ}

# 1. demux barcode R1
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Demux barcodes of R1 -> LTR only"
fastq-multx -m 1 -B ${BARCODE_LTR} ${TMPDIR}/r1.no12.fastq.gz ${TMPDIR}/r2.no12.fastq.gz -o ${TMPDIR}/r1.no12.%.fastq.gz -o ${TMPDIR}/r2.no12.%.fastq.gz

# 2. demux barcode R2 (revert pairs!!)
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Demux barcodex of R2 -> LC only"
for k in ${BCLTR[@]}; do
  fastq-multx -m 1 -B ${BARCODE_LC} ${TMPDIR}/r2.no12.${k}.fastq.gz ${TMPDIR}/r1.no12.${k}.fastq.gz -o ${TMPDIR}/bcmuxall/r2.no12.${k}.%.fastq.gz -o ${TMPDIR}/bcmuxall/r1.no12.${k}.%.fastq.gz
 # rm ${TMPDIR}/r1.no12.${k}.fastq.gz ${TMPDIR}/r2.no12.${k}.fastq.gz ${TMPDIR}/bcmuxall/r1.no12.${k}.unmatched.fastq.gz
done
for k in ${BCLTR[@]}; do
  fastq-multx -m 1 -B ${BARCODE_LC} ${TMPDIR}/r2.no12.${k}.fastq.gz ${TMPDIR}/r1.no12.${k}.fastq.gz -o ${TMPDIR}/bcmuxall/r2.no12.${k}.%.fastq.gz -o ${TMPDIR}/bcmuxall/r1.no12.${k}.%.fastq.gz
  rm ${TMPDIR}/r1.no12.${k}.fastq.gz ${TMPDIR}/r2.no12.${k}.fastq.gz ${TMPDIR}/bcmuxall/r1.no12.${k}.unmatched.fastq.gz
done
#rm ${TMPDIR}/r1.no12.fastq.gz ${TMPDIR}/r2.no12.fastq.gz ${TMPDIR}/r1.no12.unmatched.fastq.gz ${TMPDIR}/r2.no12.unmatched.fastq.gz

for b in ${BCLTR[@]}; do
  for k in ${BCLC[@]}; do
    TAG="${b}.${k}"
    # echo $TAG
    for a in ${ASSOBCLIST[@]}; do
      if [ ${TAG} == ${a} ]; then echo "--- This tag exists in the association file -> analyze it ---";
   #  		##### ================ TRIMMING ======================== #####
     		echo $TAG
        echo ${LC_rev}
			#3) trim LC from R1. If trimming this read will remove the whole read, DO NOT rescue it! (policy: everything that is trimmed in R1 BUT will lead to a zero or too small read is removed)
			#3.1) from the paired-end reads
			echo "trim lc from r1 from the paired-end reads"
			#flexbar2.5 --reads ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.fastq.gz --reads2 ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.fastq.gz --target ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LC -f i1.8 -a "${LC_rev}" --threads ${MAXTHREADS} -ae RIGHT -at 4 -ao 8 -ai -4 -m 2 -q 1 --max-uncalled 800 -z GZ
      flexbar --reads ${TMPDIR}/bcmuxall/r1.no12.${b}.${k}.fastq.gz --reads2 ${TMPDIR}/bcmuxall/r2.no12.${b}.${k}.fastq.gz --target ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LC -qf i1.8 -a ${LC_rev} --threads ${MAXTHREADS} --adapter-trim-end RIGHT -ae 0.4 -ao 8 -ai -4 -m 2 --max-uncalled 800 -z GZ
      echo "trimming 2"

			flexbar -r ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LC_2.fastq.gz -t ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC -qf i1.8 -n ${MAXTHREADS} -a ${LC_fwd} -at LEFT -ao 14 -z GZ -o --adapter-error-rate 0.1 --adapter-trimmed-out ONLY -u 8000000
			zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC.fastq.gz | fastq_to_fasta  -Q 33 -n | fasta2csv | cut -d' ' -f1 > ${TMPDIR}/bcmuxall/r2.${b}.${k}.cleanR1LCR2LC.list
			sort --parallel=5 ${TMPDIR}/bcmuxall/r2.${b}.${k}.cleanR1LCR2LC.list > ${TMPDIR}/bcmuxall/r2.${b}.${k}.cleanR1LCR2LC.sorted.list
			zcat ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LC_1.fastq.gz | python3 ${RAAVIOLIDIR}/utils_scripts/fqextract_pureheader.v3.py <(cat ${TMPDIR}/bcmuxall/r2.${b}.${k}.cleanR1LCR2LC.sorted.list) | pigz --best -f -c > ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_1.fastq.gz
			mv ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC.fastq.gz ${TMPDIR}/bcmuxall/${b}.${k}.cleanR1LCR2LC_2.fastq.gz

		  fi
	done
  done
done