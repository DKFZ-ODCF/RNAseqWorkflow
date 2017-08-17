#!/usr/bin/env bash

set -vx

source $CONFIG_FILE
source $TOOL_NAV_LIB

if [[ ${singleCellSequencingSystem} == "wafergen" || ${singleCellSequencingSystem} == "fluidigm" ]]; then
    # TODO (Jeongbin): Ultimately, it would be better to merge 'parse_barcode.py' into Roddy workflow to avoid code duplication
    #export READS_KALLISTO=(e `${PYTHON_BINARY} ${TOOL_PARSE_BARCODE_FILE} -s${SAMPLE} ${BARCODE_FILE} | cut -f 4 | awk -v P=${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}/ '{print P $0}'` )
    export READS_STAR_LEFT=`${PYTHON_BINARY} ${TOOL_PARSE_BARCODE_FILE} -s${SAMPLE} ${BARCODE_FILE} | cut -f 4 | awk -v P=${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}/ '{print P $0}' | paste -d, -s`
    export PARM_READGROUPS=`${PYTHON_BINARY} ${TOOL_PARSE_BARCODE_FILE} -t1 -s${SAMPLE} ${BARCODE_FILE}`
    bash ${TOOL_STAR_ALIGNMENT}
#elif [[ ${singleCellSequencingSystem} == "dropseq" ]]; then
#    # TODO: support dropseq
else
    echo "Not supported platform: ${singleCellSequencingSystem}"
    exit 1
fi
