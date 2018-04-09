#!/usr/bin/env bash


module load "python/${PYTHON_VERSION:?No PYTHON_VERSION}"
export PYTHON_BINARY=python

module load "STAR/${STAR_VERSION:?No STAR_VERSION}"
export STAR_BINARY=STAR

module load "subread/${SUBREAD_VERSION:?No SUBREAD_VERSION}"
export FEATURECOUNTS_BINARY=featureCounts

module load "sambamba/${SAMBAMBA_VERSION:?No SAMBAMBA_VERSION}"
export SAMBAMBA_BINARY=sambamba

module load "samtools/${SAMTOOLS_VERSION:?No SAMTOOLS_VERSION}"
export SAMTOOLS_BINARY=samtools

if [ $SAMTOOLS_VERSION = "1.3.1" ]
then
	module load "htslib/${HTSLIB_VERSION:?No HTSLIB_VERSION}"
fi

module load "rnaseqc/${RNASEQC_VERSION:?No RNASEQC_VERSION}"
export RNASEQC_BINARY=rnaseqc

module load "kallisto/${KALLISTO_VERSION:?No KALLISTO_VERSION}"
export KALLISTO_BINARY=kallisto

module load "qualimap/${QUALIMAP_VERSION:?No QUALIMAP_VERSION}"
export QUALIMAP_BINARY=qualimap

module load "arriba/${ARRIBA_VERSION:?No ARRIBA_VERSION}"
export ARRIBA_BINARY=arriba
export ARRIBA_READTHROUGH_BINARY=extract_read-through_fusions
export ARRIBA_DRAW_FUSIONS=draw_fusions.R

module load "R/${R_VERSION:?No R_VERSION}"

export JOB_PROFILER_BINARY="strace.sh"

