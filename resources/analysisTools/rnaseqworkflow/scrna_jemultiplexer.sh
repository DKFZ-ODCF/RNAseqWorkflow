#!/usr/bin/env bash

set -vx

# source $CONFIG_FILE
# source $CONFIG_FILE
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

READ_LEFT_UMI_PROCESSED="${RODDY_SCRATCH}/left_umi.fastq"
READ_RIGHT_UMI_PROCESSED="${RODDY_SCRATCH}/right_umi.fastq"

mkfifo ${READ_LEFT_COMBINED}
READS_LEFT=`echo ${READS_LEFT} | tr ',' ' '`
eval "zcat ${READS_LEFT} > ${READ_LEFT_COMBINED} &"

mkfifo ${READ_RIGHT_COMBINED}
READS_RIGHT=`echo ${READS_RIGHT} | tr ',' ' '`
eval "zcat ${READS_RIGHT} > ${READ_RIGHT_COMBINED} &"

if [[ ${processUMI} == true ]]; then
    mkfifo ${READ_LEFT_UMI_PROCESSED}
    mkfifo ${READ_RIGHT_UMI_PROCESSED}
    ${PYPY_BINARY} ${TOOL_EXTRACT_UMI} ${READ_LEFT_COMBINED} ${READ_RIGHT_COMBINED} ${READ_LEFT_UMI_PROCESSED} ${READ_RIGHT_UMI_PROCESSED} ${UMI_PATTERN} &
    READ_LEFT_JE_INPUT=${READ_LEFT_UMI_PROCESSED}
    READ_RIGHT_JE_INPUT=${READ_RIGHT_UMI_PROCESSED}
else
    READ_LEFT_JE_INPUT=${READ_LEFT_COMBINED}
    READ_RIGHT_JE_INPUT=${READ_RIGHT_COMBINED}
fi

make_directory $JE_OUTDIR
# echo_run "$JE_BINARY demultiplex FASTQ_FILE1=${READ_LEFT_JE_INPUT} FASTQ_FILE2=${READ_RIGHT_JE_INPUT} OUTPUT_DIR=${JE_OUTDIR} $JE_PARAMS"
echo_run "$JE_BINARY FASTQ_FILE1=${READ_LEFT_JE_INPUT} FASTQ_FILE2=${READ_RIGHT_JE_INPUT} OUTPUT_DIR=${JE_OUTDIR} $JE_PARAMS"

rm -f ${READ_LEFT_COMBINED}
rm -f ${READ_RIGHT_COMBINED}
if [[ ${processUMI} == true ]]; then
    rm -f ${READ_LEFT_UMI_PROCESSED}
    rm -f ${READ_RIGHT_UMI_PROCESSED}
fi

touch ${CHECKPOINT_DEMULTIPLEXING}
