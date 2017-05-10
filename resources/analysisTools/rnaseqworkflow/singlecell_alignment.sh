#!/usr/bin/env bash

set -vx

source $CONFIG_FILE
source $TOOL_NAV_LIB

# Below values will be redefined
export SAMPLE
export READS_STAR_LEFT
export READS_STAR_RIGHT
export PARM_READGROUPS
export READS_KALLISTO
export CHECKPOINT_ALIGNMENT

READS_STAR_RIGHT=
if [[ ${singleCellSequencingSystem} == "wafergen" || ${singleCellSequencingSystem} == "fluidigm" ]]; then
    SAMPLES=`${PYTHON_BINARY} ${TOOL_PARSE_BARCODE_FILE} -t2 ${BARCODE_FILE}`
    CHECKPOINT_ALIGNMENT_MERGED=${CHECKPOINT_ALIGNMENT}
    rm -f ${CHECKPOINT_ALIGNMENT_MERGED}
    for s in ${SAMPLES}; do
        SAMPLE=${s}
        READS_KALLISTO=`${PYTHON_BINARY} ${TOOL_PARSE_BARCODE_FILE} -s${SAMPLE} ${BARCODE_FILE} | cut -f 4 | awk -v P=${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}/ '{print P $0}'`
        READS_STAR_LEFT=`echo ${READS_KALLISTO} | tr ' ' ','`
        #READS_STAR_LEFT=${READS_STAR_LEFT%,} # Remove trailing comma
        PARM_READGROUPS=`${PYTHON_BINARY} ${TOOL_PARSE_BARCODE_FILE} -t1 -s${SAMPLE} ${BARCODE_FILE}`
        CHECKPOINT_ALIGNMENT=${outputAnalysisBaseDirectory}/.checkpoint_${SAMPLE}_alignment

        source $CONFIG_FILE
        bash ${TOOL_STAR_ALIGNMENT}

        cat ${CHECKPOINT_ALIGNMENT} >> ${CHECKPOINT_ALIGNMENT_MERGED}
    done
#elif [[ ${singleCellSequencingSystem} == "dropseq" ]]; then
#    # TODO: support dropseq
else
    echo "Not supported platform: ${singleCellSequencingSystem}"
    exit 1
fi
