#!/usr/bin/env bash

set -vx

source $CONFIG_FILE

if [ "$LOAD_MODULE" == true ]
then
	module load $MODULE_ENV
	CUTADAPT_BINARY=cutadapt
fi

mkdir -p ${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}
mkdir -p ${SCRATCH}/${PID}_TRIMMING

trimOpts="-b ${CUTADAPT_BINARY} -t ${SCRATCH}/${PID}_TRIMMING -a ${ADAPTER_SEQ} -c ${CORES} -q ${TRIM_3P_QUAL}"

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
        READ1_FILE=${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}/read1.txt
        ${PYTHON_BINARY} ${TOOL_PARSE_BARCODE_FILE} ${BARCODE_FILE} | cut -f 4 | awk -v P=${outputAnalysisBaseDirectory}/${demultiplexOutputDirectory}/ '{print P $0}' > ${READ1_FILE}
    #elif [[ ${singleCellSequencingSystem} == "dropseq" ]]; then
    #    # TODO: support dropseq
    fi
    trimOpts="${trimOpts} -1 ${READ1_FILE}"
else
    READ1_FILE=${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}/read1.txt
    READ2_FILE=${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}/read2.txt

    echo ${READ1} | tr ' ' '\n' > ${READ1_FILE}
    echo ${READ2} | tr ' ' '\n' > ${READ2_FILE}

    trimOpts="${trimOpts} -1 ${READ1_FILE} -2 ${READ2_FILE}"
fi

echo "Cutadapt wrapper options: ${trimOpts}"
eval "python ${TOOL_CUTADAPT_WRAPPER} ${trimOpts} ${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}"

rm -f $READ1_FILE
if [ ! -z $READ2_FILE ]; then
    rm -f $READ2_FILE
fi
rm -rf ${SCRATCH}/${PID}_TRIMMING

touch ${CHECKPOINT_TRIMMING}