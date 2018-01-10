#!/usr/bin/env bash

set -vx

# source $CONFIG_FILE
# source $CONFIG_FILE

if [ "$RUN_FEATURE_COUNTS" == true ]
then
    cd $COUNT_DIR
    if [ `ls ${SAMPLE}_${PID}_*.s0 | wc -l` == $BAM_FILE_NUM ]; then
        FC_OUTS=
        for (( i=0; i<$BAM_FILE_NUM; i++ )); do
            FC_OUTS="${FC_OUTS} ${SAMPLE}_${PID}_${i}.featureCounts.s0"
        done
        ${PYTHON_BINARY} ${TOOL_MERGE_FEATURECOUNTS_TABLES} ${GENE_MODELS} ${FC_OUTS}
    else
        echo "Error: not all featureCounts output files are available."
        exit 1
    fi
fi

if [ "$RUN_LIBRARY_QC" == true ]
then
    mkdir -p ${LIBRARY_QC_DIR}
    if [ `ls ${SAMPLE}_${PID}_*.s0.summary | wc -l` == $BAM_FILE_NUM ]; then
        FC_SUMMARIES=
        for (( i=0; i<$BAM_FILE_NUM; i++ )); do
            FC_SUMMARIES="${FC_SUMMARIES} ${SAMPLE}_${PID}_${i}.featureCounts.s0.summary"
        done
        ${PYTHON_BINARY} ${TOOL_SINGLE_CELL_LIBRARY_QC} ${SAMPLE} ${PID} ${COUNT_DIR} ${LIBRARY_QC_DIR} ${FC_SUMMARIES}
    else
        echo "Error: not all featureCounts summary files are available."
        exit 2
    fi
fi

touch ${CHECKPOINT_POSTPROCESSING}