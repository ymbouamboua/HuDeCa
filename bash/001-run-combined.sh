#!/bin/bash

#SAMPLE="FN_S1256"
SAMPLE="FN_S3478"
WD="/data/analysis/data_mbouamboua"
DEMUXAFY="/data/analysis/data_mbouamboua/Demuxafy"
OUTDIR=${WD}/data_analysis/hudeca/WS_demuxafy/${SAMPLE}
OUTDIR_MAJORITY=${WD}/data_analysis/hudeca/WS_demuxafy/${SAMPLE}/combine_majoritySinglet
OUTDIR_ANYDOUBLET=${WD}/data_analysis/hudeca/WS_demuxafy/${SAMPLE}/combine_anyDoublet

if [ ! -d ${OUTDIR_MAJORITY} ]; then
	mkdir ${OUTDIR_MAJORITY}
    echo "$OUTDIR_MAJORITY created" ;
fi

if [ ! -d ${OUTDIR_ANYDOUBLET} ]; then
	mkdir ${OUTDIR_ANYDOUBLET}
    echo "$OUTDIR_ANYDOUBLET created" ;
fi

METHOD="MajoritySinglet"
singularity exec --bind ${WD} ${DEMUXAFY}/Demuxafy.sif Combine_Results.R -o $OUTDIR_MAJORITY/combined_results.tsv --vireo ${OUTDIR}/Vireo_OUTDIR --souporcell ${OUTDIR}/Souporcell_OUTDIR --scds ${OUTDIR}/Scds_OUTDIR --method ${METHOD}

METHOD="AnyDoublet"
singularity exec --bind ${WD} ${DEMUXAFY}/Demuxafy.sif Combine_Results.R -o $OUTDIR_ANYDOUBLET/combined_results.tsv --vireo ${OUTDIR}/Vireo_OUTDIR --souporcell ${OUTDIR}/Souporcell_OUTDIR --scds ${OUTDIR}/Scds_OUTDIR --method ${METHOD}

# singularity exec --bind ${WD} ${DEMUXAFY}/Demuxafy.sif Combine_Results.R -o $OUTDIR_assign/combined_results.tsv --vireo ${OUTDIR}/Vireo_OUTDIR --souporcell_assignments ${OUTDIR}/Souporcell_OUTDIR --method ${METHOD}
