#!/usr/bin/env bash

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
    export ARRIBA_DRAW_FUSIONS=draw_fusions.R
fi
