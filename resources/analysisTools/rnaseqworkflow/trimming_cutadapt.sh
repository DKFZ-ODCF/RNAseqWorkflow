#!/usr/bin/env bash

set -vx

source $CONFIG_FILE

if [ "$LOAD_MODULE" == true ]
then
	module load $MODULE_ENV
	CUTADAPT_BINARY=cutadapt
fi

mkdir -p ${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}

if [ -z "${CHUNK_INDEX}" ]; then
    CUTADAPTTMPDIR=${SCRATCH}/${PID}_${SAMPLE}_TRIMMING
else
    CUTADAPTTMPDIR=${SCRATCH}/${PID}_${SAMPLE}_${CHUNK_INDEX}_TRIMMING
fi
mkdir -p ${CUTADAPTTMPDIR}

trimOpts="-b ${CUTADAPT_BINARY} -t ${CUTADAPTTMPDIR} -a ${ADAPTER_SEQ} -c ${CORES} -q ${TRIM_3P_QUAL}"

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

READ1_FILE=${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}/read1.txt
echo ${READ1} | tr ' ' '\n' > ${READ1_FILE}

if [ "$runSingleCellWorkflow" == true ] || [ "$useSingleEndProcessing" == true ]; then
    trimOpts="${trimOpts} -1 ${READ1_FILE}"
else
    READ2_FILE=${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}/read2.txt
    echo ${READ2} | tr ' ' '\n' > ${READ2_FILE}
    trimOpts="${trimOpts} -1 ${READ1_FILE} -2 ${READ2_FILE}"
fi

echo "Cutadapt wrapper options: ${trimOpts}"
eval "python ${TOOL_CUTADAPT_WRAPPER} ${trimOpts} ${outputAnalysisBaseDirectory}/${trimmingOutputDirectory}"

rm -f $READ1_FILE
if [ ! -z $READ2_FILE ]; then
    rm -f $READ2_FILE
fi
rm -rf ${CUTADAPTTMPDIR}

touch ${CHECKPOINT_TRIMMING}