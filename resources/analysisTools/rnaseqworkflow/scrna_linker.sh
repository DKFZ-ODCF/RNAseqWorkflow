#!/usr/bin/env bash

set -vx

OUTDIR=${outputAnalysisBaseDirectory}/${demultiplexOutputDirectory}

mkdir -p $OUTDIR
cd $OUTDIR

for tgt in ${SYMLINK_TARGETS}; do
    LNCMD="ln -s `echo ${tgt} | tr ',' ' '`"
    eval $LNCMD
done

touch ${CHECKPOINT_LINKING}
