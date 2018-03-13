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
    JE_PARAMS="${JE_PARAMS_WAFERGEN} OUTPUT_DIR=${JE_OUTDIR}"
#elif [[ ${singleCellSequencingSystem} == "dropseq" ]]; then
#    JE_PARAMS="" # TODO: set appropriate demultiplexing params for dropseq
else
    echo "Not supported platform: ${singleCellSequencingSystem}"
    exit 1
fi

nlines=`wc -l < ${BARCODE_JE}`
chunksize=100
nchunks=$(( (${nlines}+${chunksize}-1)/${chunksize} ))

if [[ -f ${JE_OUTDIR}/resume.txt ]]; then
    startidx=`cat ${JE_OUTDIR}/resume.txt`
    rm -f ${JE_OUTDIR}/unassigned_chunk$(($startidx+1))_1.fasta.gz ${JE_OUTDIR}/unassigned_chunk$(($startidx+1))_2.fasta.gz
    CLEANUP=false
else
    startidx=1

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
    CLEANUP=true

    make_directory $JE_OUTDIR
    BARCODE_FILE=${JE_OUTDIR}/jemultiplexer_chunk1.in
    head -n ${chunksize} ${BARCODE_JE} > ${BARCODE_FILE}
    echo_run "$JE_BINARY BARCODE_FILE=${BARCODE_FILE} FASTQ_FILE1=${READ_LEFT_JE_INPUT} FASTQ_FILE2=${READ_RIGHT_JE_INPUT} UNASSIGNED_FILE_NAME_1=unassigned_chunk1_1.fasta UNASSIGNED_FILE_NAME_2=unassigned_chunk1_2.fasta METRICS_FILE_NAME=jemultiplexer_out_chunk1_stats.txt $JE_PARAMS"
fi

for ((i=startidx;i<nchunks;i++)); do
    BARCODE_FILE=${JE_OUTDIR}/jemultiplexer_chunk$(($i+1)).in
    tail -n+$((${i}*${chunksize}+1)) ${BARCODE_JE} | head -n ${chunksize} > ${BARCODE_FILE}
    echo_run "$JE_BINARY BARCODE_FILE=${BARCODE_FILE} FASTQ_FILE1=${JE_OUTDIR}/unassigned_chunk${i}_1.fasta.gz FASTQ_FILE2=${JE_OUTDIR}/unassigned_chunk${i}_2.fasta.gz UNASSIGNED_FILE_NAME_1=unassigned_chunk$(($i+1))_1.fasta UNASSIGNED_FILE_NAME_2=unassigned_chunk$(($i+1))_2.fasta METRICS_FILE_NAME=jemultiplexer_out_chunk$(($i+1))_stats.txt $JE_PARAMS"
    rm -f ${JE_OUTDIR}/unassigned_chunk${i}_1.fasta.gz ${JE_OUTDIR}/unassigned_chunk${i}_2.fasta.gz
    echo ${i} > ${JE_OUTDIR}/resume.txt
done

mv ${JE_OUTDIR}/unassigned_chunk${nchunks}_1.fasta.gz ${JE_OUTDIR}/unassigned_1.fasta.gz
mv ${JE_OUTDIR}/unassigned_chunk${nchunks}_2.fasta.gz ${JE_OUTDIR}/unassigned_2.fasta.gz

if [[ "$CLEANUP" == true ]]; then
    rm -f ${READ_LEFT_COMBINED}
    rm -f ${READ_RIGHT_COMBINED}
    if [[ "$processUMI" == true ]]; then
        rm -f ${READ_LEFT_UMI_PROCESSED}
        rm -f ${READ_RIGHT_UMI_PROCESSED}
    fi
fi

touch ${CHECKPOINT_DEMULTIPLEXING}

exit 0
