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

    # rest of vispa2 block unchanged …
    # (exons_true / else branches remain identical)

    echo "

    ---------------------------------------------------------------------------------
                          ENDING PROCESSING AT: `date +'%Y-%m-%d %H:%M:%S'`
    ---------------------------------------------------------------------------------
        "
fi

if [ ${TYPE} = "vispa" ]; then
    printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Preparing Data.\n"
    case "${MATRIX}" in
    *.gz )
        N_MATRIX="`basename ${MATRIX}`";
        zcat ${MATRIX} | sed '/^M/ d' > ${N_MATRIX:0:-7}.noChrM.tsv
        zcat ${N_MATRIX:0:-7}.noChrM.tsv | sed "s/\t0/\t/g" > ${N_MATRIX:0:-7}.no0.tsv;
        rm ${N_MATRIX:0:-7}.noChrM.tsv;
        MATRIX="`basename ${N_MATRIX:0:-7}.no0.tsv`";
        ;;
    *)
        N_MATRIX="`basename ${MATRIX}`";
        cat ${MATRIX} | sed '/^M/ d' > ${N_MATRIX:0:-4}.noChrM.tsv
        sed "s/\t0/\t/g" ${N_MATRIX:0:-4}.noChrM.tsv > ${N_MATRIX:0:-4}.no0.tsv;
        MATRIX="`basename ${N_MATRIX:0:-4}.no0.tsv`";
        ;;
    esac

    printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Annotation.\n"
   python3 ${RAAVIOLIDIR}/utils_scripts/annotate_bed.v3.py  -a ${GTF} -b ${MATRIX} --matrix -o ${OUTDIR}/iss.annotation.tsv

    # rest of vispa block unchanged …
    # (exons_true / else branches remain identical)

    echo "

    ---------------------------------------------------------------------------------
                          ENDING PROCESSING AT: `date +'%Y-%m-%d %H:%M:%S'`
    ---------------------------------------------------------------------------------
        "
fi

if [ ${TYPE} = "clustering" ]; then
    printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Preparing Data.\n"
    case "${MATRIX}" in
    *.gz )
        N_MATRIX="`basename ${MATRIX}`";
        zcat ${MATRIX} | sed 's/,$//g' > ${N_MATRIX:0:-7}.tmp.csv;
        sed 's|,|\t|g' ${N_MATRIX:0:-7}.tmp.csv > ${N_MATRIX:0:-7}.converted.tsv;
        rm ${N_MATRIX:0:-7}.tmp.csv;
        sed "s/\t0/\t/g" ${N_MATRIX:0:-7}.converted.tsv > ${N_MATRIX:0:-7}.no0.tsv;
        rm ${N_MATRIX:0:-7}.converted.tsv;
        MATRIX="`basename ${N_MATRIX:0:-7}.no0.tsv`";
        cat ${MATRIX} | awk -F$'\t' '{print $2"\t"$3"\t"$3}' | sed '1d' > ${MATRIX:0:-7}.bed

        printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Annotation.\n"
        python3 ${RAAVIOLIDIR}/utils_scripts/annotate_bed.v3.py -a ${GTF} -b ${MATRIX:0:-7}.bed -o ${OUTDIR}/iss.annotation.tsv
        rm ${MATRIX:0:-7}.bed;
        ;;
    *)
        N_MATRIX="`basename ${MATRIX}`";
        sed 's/,$//g' ${MATRIX} > ${N_MATRIX:0:-4}.tmp.csv;
        sed 's|,|\t|g' ${N_MATRIX:0:-4}.tmp.csv > ${N_MATRIX:0:-4}.converted.tsv;
        rm ${N_MATRIX:0:-4}.tmp.csv;
        sed "s/\t0/\t/g" ${N_MATRIX:0:-4}.converted.tsv > ${N_MATRIX:0:-4}.no0.tsv;
        rm ${N_MATRIX:0:-4}.converted.tsv;
        MATRIX="`basename ${N_MATRIX:0:-4}.no0.tsv`";
        cat ${MATRIX} | awk -F$'\t' '{print $2"\t"$3"\t"$3}' | sed '1d' > ${MATRIX:0:-4}.bed

        printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Annotation.\n"
        python3 ${RAAVIOLIDIR}/utils_scripts/annotate_bed.v3.py -a ${GTF} -b ${MATRIX:0:-4}.bed -o ${OUTDIR}/iss.annotation.tsv
        rm ${MATRIX:0:-4}.bed;
        ;;
    esac

    # rest of clustering block unchanged …
    # (exons_true / else branches remain identical)

    echo "

    ---------------------------------------------------------------------------------
                          ENDING PROCESSING AT: `date +'%Y-%m-%d %H:%M:%S'`
    ---------------------------------------------------------------------------------
        "
fi
