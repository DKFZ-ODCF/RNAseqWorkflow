#!/usr/bin/env bash

module load ""
export STAR_BINARY=STAR
export FEATURECOUNTS_BINARY=featureCounts
export SAMBAMBA_BINARY=sambamba
export SAMTOOLS_BINARY=samtools # samtools_bin defaults to v0.0.19 despite the RNAseq.XML!!!!!!!!!!!!!!!!!!!!!!!!!!
export RNASEQC_BINARY=rnaseqc
export KALLISTO_BINARY=kallisto
export QUALIMAP_BINARY=qualimap
export ARRIBA_BINARY=arriba
export ARRIBA_READTHROUGH_BINARY=extract_read-through_fusions
export ARRIBA_DRAW_FUSIONS=draw_fusions.R

