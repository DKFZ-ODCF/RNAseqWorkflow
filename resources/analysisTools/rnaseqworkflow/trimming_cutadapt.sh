#!/usr/bin/env bash

set -vx

source $CONFIG_FILE

if [ "$LOAD_MODULE" == true ]
then
	module load $MODULE_ENV
	CUTADAPT_BINARY=cutadapt
fi

mkdir -p ${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}

trimOpts="-b ${CUTADAPT_BINARY} -a ${ADAPTER_SEQ} -c ${CORES} -q ${TRIM_3P_QUAL}"

if [ "$TRIM_NEXTSEQ" == true ]; then
    trimOpts="${trimOpts} -n"
fi

if [ "$TRIM_3P_POLYA" == true ]; then
    trimOpts="${trimOpts} -A3"
fi

if [ "$TRIM_3P_POLYT" == true ]; then
    trimOpts="${trimOpts} -T3"
fi

if [ "$TRIM_3P_POLYG" == true ]; then
    trimOpts="${trimOpts} -G3"
fi

if [ "$TRIM_3P_POLYN" == true ]; then
    trimOpts="${trimOpts} -N3"
fi

if [ "$runSingleCellWorkflow" == true ]; then
    if [[ ${singleCellSequencingSystem} == "wafergen" || ${singleCellSequencingSystem} == "fluidigm" ]]; then
        READ1=`${PYTHON_BINARY} ${TOOL_PARSE_BARCODE_FILE} ${BARCODE_FILE} | cut -f 4 | awk -v P=${outputAnalysisBaseDirectory}/${demultiplexOutputDirectory}/ '{print P $0}'`
        READ1=`echo ${READ1} | tr ' ' ','`
    #elif [[ ${singleCellSequencingSystem} == "dropseq" ]]; then
    #    # TODO: support dropseq
    fi
    trimOpts="${trimOpts} -1 ${READ1}"
else
    trimOpts="${trimOpts} -1 ${READ1} -2 ${READ2}"
fi

echo "Cutadapt options: ${trimOpts}"
eval "python ${TOOL_CUTADAPT_WRAPPER} ${trimOpts} ${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}"

touch ${CHECKPOINT_TRIMMING}