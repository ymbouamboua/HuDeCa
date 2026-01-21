#!/bin/bash
# Should be run on bego.ipmc.cnrs.fr
# sh Run_Demuxafy_HuDeCa.sh ../data_analysis/hudeca/data/FN_S1256/FN_S1256.txt


# Job variables
WD="/data/analysis/data_mbouamboua"
DEMUXAFY="/data/analysis/data_mbouamboua/Demuxafy"
INDIR=${WD}/data_analysis/hudeca/data
OUTDIR=${WD}/data_analysis/hudeca
OUT=${OUTDIR}/WS_demuxafy
TMP=${OUTDIR}/tmp
THREADS=30
# exp_design=$1

# Global variables
VCF=${WD}/data_analysis/genome/GRCh38/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf
VCF=${WD}/data_analysis/genome/bulk_rna.vcf.gz
GENOME=${WD}/data_analysis/genome/GRCh38/hg38.fa

# Create necessary directories if they do not exist
if [ ! -d ${TMP} ]; then
    mkdir -p ${TMP}
    echo "$TMP created" ;
fi

if [ ! -d ${OUT} ]; then
    mkdir -p ${OUT}
    echo "$OUT created" ;
fi

# Read the list of samples and process each line
while IFS= read -r line; do
    # Extract sample and tool names from the line
    SAMPLE=$(echo "$line" | awk '{print $1}')
    TOOL=$(echo "$line" | awk '{print $2}')

    # Prepare directories for each sample and tool
    OUT_SAMPLE=${OUT}/${SAMPLE}
    if [ ! -d ${OUT_SAMPLE} ]; then
        mkdir -p ${OUT_SAMPLE}
        echo "$OUT_SAMPLE created" ;
    fi
    
    OUT_TOOL=${OUT_SAMPLE}/${TOOL}_OUTDIR
    if [ ! -d ${OUT_TOOL} ]; then
        mkdir -p ${OUT_TOOL}
        echo "$OUT_TOOL created" ;
    fi

    # Define data variables for the sample
    SAMPLE_DIR=${INDIR}/${SAMPLE}
    BAM=${SAMPLE_DIR}/${SAMPLE}.bam
    BARCODES=${SAMPLE_DIR}/${SAMPLE}.barcodes
    BC_MATRIX=${SAMPLE_DIR}/${SAMPLE}_bc_matrix/
    NB_SAMPLES=$(cat ${SAMPLE_DIR}/${SAMPLE}.nbsamples | tr -d '\n')

    # Prepare and submit the job
    JOB_NAME=${SAMPLE}_${TOOL}
    LOG=${TMP}/${JOB_NAME}.out
    SCRIPT=${DEMUXAFY}/TOOLS/Run_${TOOL}.sh
    COMMAND="$SCRIPT $WD $DEMUXAFY $BAM $BARCODES $GENOME $OUT_TOOL $NB_SAMPLES $VCF $THREADS $BC_MATRIX"

    echo "# Submit ${JOB_NAME}:"
    echo "Run nohup -> ${COMMAND}"
    nohup sh ${COMMAND} > $LOG &
    
done < "$1"
