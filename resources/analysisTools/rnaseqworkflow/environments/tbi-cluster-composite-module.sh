#!/usr/bin/env bash

module load "python/${PYTHON_VERSION:?No PYTHON_VERSION}"
if [ -f ${PYTHON2_VENV:?PYTHON2_VENV variable is undefined} ]; then
    source "${PYTHON2_VENV}"
fi
export PYTHON_BINARY=python
export CUTADAPT_BINARY=cutadapt

if [ "$LOAD_MODULE" == true ]; then
    module load "${MODULE_ENV:?MODULE_ENV variable is undefined}"
    export STAR_BINARY=STAR
    export FEATURECOUNTS_BINARY=featureCounts
    export SAMBAMBA_BINARY=sambamba
    export SAMTOOLS_BINARY=samtools # samtools_bin defaults to v0.0.19 despite the RNAseq.XML!!!!!!!!!!!!!!!!!!!!!!!!!!
    export RNASEQC_BINARY=rnaseqc
    export KALLISTO_BINARY=kallisto
    export QUALIMAP_BINARY=qualimap
    export ARRIBA_BINARY=arriba
    export RSCRIPT_BINARY=Rscript
    export ARRIBA_READTHROUGH_BINARY=extract_read-through_fusions
    export ARRIBA_DRAW_FUSIONS=draw_fusions.R
fi