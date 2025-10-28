#!/bin/bash

echo "
        +--------------------------------------------------------+
        |                                                        |
        |               ISs Matrix Annotation V2                 |
        |                                                        |
        +--------------------------------------------------------+
        |  Author:   Giulio Spinozzi, PhD                        |
        |  Date:     April 2016                              	 |
        |  Version:  2.1                                         |
        |  Contact:  spinozzi.giulio@hsr.it                      |
        +--------------------------------------------------------+
"


##### ============================== RUN INFO =============================== #####
RUN_STARTED_AT=`date +"%Y-%m-%d %H:%M:%S"`;
RUN_ID="`whoami`"" ${RUN_STARTED_AT}";

TODAY=`date +"%Y%m%d%H%M"`;
#=================================================================================#



##### ================================ ARGS ================================= #####
usage()
{
    echo "This app annotate ISs matrix file (csv, tsv)."
    echo
    echo "Usage: $0 [-m ISs.matrix.tsv] [-t type] [-g gtf_file.gtf] [-o output dir]"
    echo
    echo "  [-m ISs.matrix.tsv]   -   Input File: ISs Matrix (CSV, TSV)"
    echo "  [-t type]             -   Input Matrix Type: vispa2, vispa or clustering"
    echo "  [-e exons]            -   Exon Info: exons_true or exons_false. Default False"
    echo "  [-g gtf_file.gtf]     -   Input File: GTF File"
    echo "  [-o output dir]       -   Output Directory"
    echo
    exit
}

while getopts ":g:m:o:t:e:h" Option
    do
    case $Option in
        g ) GTF=$OPTARG ;;
        m ) MATRIX=$OPTARG ;;
        o ) OUTDIR=$OPTARG ;;
        t ) TYPE=$OPTARG ;;
        e ) EXONS=$OPTARG ;;
        h ) usage ;;
        * ) echo "unrecognized argument. use '-h' for usage information."; exit -1 ;;
    esac
done
shift $(($OPTIND - 1))
#=================================================================================#



##### ========================== PRELIMINARY CHECKS ========================= #####
## Arguments
if [ -z "$MATRIX" ]; then
    usage
fi

if [ -z "$GTF" ]; then
    usage
fi

if [ -z "$OUTDIR" ]; then
    usage
fi

if [ -z "$TYPE" ]; then
    usage
fi

if [ ! -f ${MATRIX} ]; then
    echo ""
    echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    echo " ${MATRIX} NOT EXISTS!!!!! "
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    exit 0
fi

if [ ! -f ${GTF} ]; then
    echo ""
    echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    echo " ${GTF} NOT EXISTS!!!!! "
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    exit 0
fi
#=================================================================================#


#---------------------------------------***---------------------------------------#
##### =============================== PROGRAM =============================== #####
echo "

---------------------------------------------------------------------------------
                    STARTING PROCESSING AT: $RUN_STARTED_AT
---------------------------------------------------------------------------------
    "


##### =========================== SETUP PARAMETERS ========================== #####
mkdir ${OUTDIR};
#=================================================================================#

if [ ${TYPE} = "vispa2" ]; then
    printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Preparing Data.\n"
    printf "chr\tintegration_locus\tstrand\n" > IS.positions.tsv;
    case "${MATRIX}" in
    *.gz )
        N_MATRIX="`basename ${MATRIX}`";
        zcat ${MATRIX} | sed '/^chrM/ d' > ${N_MATRIX:0:-7}.noChrM.tsv
        cut -f1 ${N_MATRIX:0:-7}.noChrM.tsv | tail -n +2 | sed "s/\_/\t/g" | sed 's|chr||g' >> IS.positions.tsv;
        cut -f2- ${N_MATRIX:0:-7}.noChrM.tsv > IS.quantification.tsv
        if [ -z $(grep -P "\t0\t" IS.quantification.tsv) ]; then
            paste -d '\t' IS.positions.tsv IS.quantification.tsv | sed "s/\t0.0/\t/g" > ${N_MATRIX:0:-7}.no0.tsv;
        else
            paste -d '\t' IS.positions.tsv IS.quantification.tsv | sed "s/\t0/\t/g"  > ${N_MATRIX:0:-7}.no0.tsv;
        fi
        rm ${N_MATRIX:0:-7}.noChrM.tsv;
        MATRIX="`basename ${N_MATRIX:0:-7}.no0.tsv`";
        ;;
    *)
        N_MATRIX="`basename ${MATRIX}`";
        cat ${MATRIX} | sed '/^chrM/ d' > ${N_MATRIX:0:-4}.noChrM.tsv
        cut -f1 ${N_MATRIX:0:-4}.noChrM.tsv | tail -n +2 | sed "s/\_/\t/g" | sed 's|chr||g' >> IS.positions.tsv;
        cut -f2- ${N_MATRIX:0:-4}.noChrM.tsv > IS.quantification.tsv;
        if [ -z $(grep -P "\t0\t" IS.quantification.tsv) ]; then
            paste -d '\t' IS.positions.tsv IS.quantification.tsv | sed "s/\t0.0/\t/g" > ${N_MATRIX:0:-4}.no0.tsv;
        else
            paste -d '\t' IS.positions.tsv IS.quantification.tsv | sed "s/\t0/\t/g"  > ${N_MATRIX:0:-4}.no0.tsv;
        fi
        rm ${N_MATRIX:0:-4}.noChrM.tsv;
        MATRIX="`basename ${N_MATRIX:0:-4}.no0.tsv`";
        ;;
    esac

    rm IS.positions.tsv IS.quantification.tsv;

    printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Annotation.\n"
    python3 ${RAAVIOLIDIR}/utils_scripts/annotate_bed.v3.py -a ${GTF} -b ${MATRIX} --matrix -o ${OUTDIR}/iss.annotation.tsv

    rm None;
    if [[ ${EXONS} = "exons_true" ]]; then
        printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Retrieving Annotation Data.\n\n"
        cat ${OUTDIR}/iss.annotation.tsv | cut -f7,11 > ${OUTDIR}/genes_data.tsv
        cat ${OUTDIR}/iss.annotation.tsv | cut -f14 > ${OUTDIR}/exons_data2.tsv
        cat ${OUTDIR}/iss.annotation.tsv | cut -f13 | awk -F '";' '{print $3}' | tr -d " " | sed 's/\bgene[^ ]*//' | awk -F '"' '{print $2}' > ${OUTDIR}/exons_data1.tsv
        cat ${OUTDIR}/iss.annotation.tsv | cut -f13 | awk -F '";' '{print $2}' | tr -d " " | awk -F '"' '{print $2}' > ${OUTDIR}/transcripts_data.tsv
        printf "GeneName\tGeneStrand\n" > ${OUTDIR}/genes.tsv
        printf "ExonNumber\tExonDistance\n" > ${OUTDIR}/exons.tsv
        printf "Transcript\n" > ${OUTDIR}/transcripts.tsv
        cat ${OUTDIR}/genes_data.tsv >> ${OUTDIR}/genes.tsv
        cat ${OUTDIR}/transcripts_data.tsv >> ${OUTDIR}/transcripts.tsv
        paste -d '\t' ${OUTDIR}/exons_data1.tsv ${OUTDIR}/exons_data2.tsv > ${OUTDIR}/exons_data.tsv
        rm ${OUTDIR}/exons_data1.tsv ${OUTDIR}/exons_data2.tsv
        cat ${OUTDIR}/exons_data.tsv >> ${OUTDIR}/exons.tsv

        printf "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Merging Annotation Data into ISs Matrix.\n"
        n_columns=$(head -n 1 ${MATRIX} | awk '{print NF}')
        cat ${MATRIX} | cut -f1,2,3 > ${OUTDIR}/matrix.part1.tsv
        cat ${MATRIX} | cut -f4-${n_columns} > ${OUTDIR}/matrix.part2.tsv
        MATRIX="`basename ${MATRIX:0:-4}`";
        paste -d '\t' ${OUTDIR}/matrix.part1.tsv ${OUTDIR}/genes.tsv ${OUTDIR}/exons.tsv ${OUTDIR}/transcripts.tsv > ${OUTDIR}/matrix.part1.annotated.tsv
        paste -d '\t' ${OUTDIR}/matrix.part1.annotated.tsv ${OUTDIR}/matrix.part2.tsv > ${OUTDIR}/${MATRIX}.annotated.tsv
        sed -i -e "s/\r//g" ${OUTDIR}/${MATRIX}.annotated.tsv
        pigz -f ${OUTDIR}/${MATRIX}.annotated.tsv
        rm ${OUTDIR}/matrix.part1.tsv ${OUTDIR}/matrix.part2.tsv ${OUTDIR}/genes.tsv ${OUTDIR}/genes_data.tsv ${OUTDIR}/exons.tsv ${OUTDIR}/exons_data.tsv ${OUTDIR}/matrix.part1.annotated.tsv ${OUTDIR}/iss.annotation.tsv ${MATRIX:0:-4}.no0.tsv ${OUTDIR}/transcripts_data.tsv ${OUTDIR}/transcripts.tsv;
    else
        printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Retrieving Annotation Data.\n\n"
        cat ${OUTDIR}/iss.annotation.tsv | cut -f7,11 > ${OUTDIR}/genes_data.tsv
        printf "GeneName\tGeneStrand\n" > ${OUTDIR}/genes.tsv
        cat ${OUTDIR}/genes_data.tsv >> ${OUTDIR}/genes.tsv

        printf "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Merging Annotation Data into ISs Matrix.\n"
        n_columns=$(head -n 1 ${MATRIX} | awk '{print NF}')
        cat ${MATRIX} | cut -f1,2,3 > ${OUTDIR}/matrix.part1.tsv
        cat ${MATRIX} | cut -f4-${n_columns} > ${OUTDIR}/matrix.part2.tsv
        MATRIX="`basename ${MATRIX:0:-4}`";
        paste -d '\t' ${OUTDIR}/matrix.part1.tsv ${OUTDIR}/genes.tsv > ${OUTDIR}/matrix.part1.annotated.tsv
        paste -d '\t' ${OUTDIR}/matrix.part1.annotated.tsv ${OUTDIR}/matrix.part2.tsv > ${OUTDIR}/${MATRIX}.annotated.tsv
        pigz -f ${OUTDIR}/${MATRIX}.annotated.tsv
        rm ${OUTDIR}/matrix.part1.tsv ${OUTDIR}/matrix.part2.tsv ${OUTDIR}/genes.tsv ${OUTDIR}/genes_data.tsv ${OUTDIR}/matrix.part1.annotated.tsv ${OUTDIR}/iss.annotation.tsv ${MATRIX:0:-4}.no0.tsv;
    fi

    echo "

    ---------------------------------------------------------------------------------
                          ENDING PROCESSING AT: `date +'%Y-%m-%d %H:%M:%S'`
    ---------------------------------------------------------------------------------
        "
fi

