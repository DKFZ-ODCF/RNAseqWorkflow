#!/usr/bin/env bash

set -vx

source $CONFIG_FILE
source $TOOL_NAV_LIB

JE_OUTDIR=${outputAnalysisBaseDirectory}/${demultiplexOutputDirectory}

echo "$READS_LEFT"
echo "$READS_RIGHT"
echo "$BARCODE_JE"
exit 0

if [ "$LOAD_MODULE" == true ]
then
	module load $MODULE_ENV
    JE_BINARY="jemultiplexer"
fi

#check_executable "$JE_BINARY"

if [[ ${singleCellSequencingSystem} == "wafergen" || ${singleCellSequencingSystem} == "fluidigm" ]]; then
    JE_PARAMS="BARCODE_FILE=${BARCODE_JE} BARCODE_READ_POS=READ_1 BARCODE_FOR_SAMPLE_MATCHING=BOTH REDUNDANT_BARCODES=true STRICT=false MAX_MISMATCHES=1 MIN_MISMATCH_DELTA=1 MIN_BASE_QUALITY=10 XTRIMLEN=0 ZTRIMLEN=0 CLIP_BARCODE=true ADD_BARCODE_TO_HEADER=false QUALITY_FORMAT=Standard KEEP_UNASSIGNED_READ=true FORCE=true    UNASSIGNED_FILE_NAME_1=unassigned_1.txt UNASSIGNED_FILE_NAME_2=unassigned_2.txt METRICS_FILE_NAME=jemultiplexer_out_stats.txt GZIP_OUTPUTS=true WRITER_FACTORY_USE_ASYNC_IO=true STATS_ONLY=false VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false"
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
