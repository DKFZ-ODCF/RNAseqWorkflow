#!/usr/bin/env bash

# Copyright (c) 2017
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
# Load a Conda environment.

source activate "${condaEnvironmentName:?No Conda environment name defined. Please set 'condaEnvironmentName'.}" \
    || (echo "Could not load Conda environment '$condaEnvironmentName'" && exit 100)

export PYTHON_BINARY=python
export STAR_BINARY=STAR
export FEATURECOUNTS_BINARY=featureCounts
export SAMBAMBA_BINARY=sambamba
export SAMTOOLS_BINARY=samtools
export RNASEQC_BINARY=rnaseqc
export KALLISTO_BINARY=kallisto
export QUALIMAP_BINARY=qualimap
export ARRIBA_BINARY=arriba
export ARRIBA_READTHROUGH_BINARY=extract_read-through_fusions
export ARRIBA_DRAW_FUSIONS=draw_fusions.R

export JOB_PROFILER_BINARY="strace.sh"

