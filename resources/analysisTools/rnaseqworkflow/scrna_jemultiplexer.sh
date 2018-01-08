#!/usr/bin/env bash

set -vx

source $TOOL_BASH_LIB

JE_OUTDIR=${outputAnalysisBaseDirectory}/${demultiplexOutputDirectory}

if [ "$LOAD_MODULE" == true ]
then
	module load $MODULE_ENV
fi

#check_executable "$JE_BINARY"

if [[ ${singleCellSequencingSystem} == "wafergen" ]]; then
    JE_PARAMS="BARCODE_FILE=${BARCODE_JE} ${JE_PARAMS_WAFERGEN}"
#elif [[ ${singleCellSequencingSystem} == "dropseq" ]]; then
#    JE_PARAMS="" # TODO: set appropriate demultiplexing params for dropseq
else
    echo "Not supported platform: ${singleCellSequencingSystem}"
    exit 1
fi

READ_LEFT_COMBINED="${RODDY_SCRATCH}/left.fastq"
READ_RIGHT_COMBINED="${RODDY_SCRATCH}/right.fastq"

mkfifo ${READ_LEFT_COMBINED}
READS_LEFT=`echo ${READS_LEFT} | tr ',' ' '`
eval "zcat ${READS_LEFT} > ${READ_LEFT_COMBINED} &"

mkfifo ${READ_RIGHT_COMBINED}
READS_RIGHT=`echo ${READS_RIGHT} | tr ',' ' '`
eval "zcat ${READS_RIGHT} > ${READ_RIGHT_COMBINED} &"

make_directory $JE_OUTDIR
echo_run "$JE_BINARY FASTQ_FILE1=${READ_LEFT_COMBINED} FASTQ_FILE2=${READ_RIGHT_COMBINED} OUTPUT_DIR=${JE_OUTDIR} $JE_PARAMS"

rm -f ${READ_LEFT_COMBINED}
rm -f ${READ_RIGHT_COMBINED}

touch ${CHECKPOINT_DEMULTIPLEXING}
