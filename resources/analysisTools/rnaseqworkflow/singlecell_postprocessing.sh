export SAMPLE

source $CONFIG_FILE

echo "$BAM_FILE_NUM"
echo "$SAMPLE"
exit 0

if [ "$RUN_FEATURE_COUNTS" == true ]
then
    cd $COUNT_DIR
    for SAMPLE in ${SAMPLES}; do
        ${PYTHON_BINARY} ${TOOL_MERGE_FEATURECOUNTS_TABLES} ${SAMPLE}_${PID}_*.featureCounts.tsv
    done
fi

if [ "$RUN_FEATURE_COUNTS_DEXSEQ" == true ]
then
    cd $COUNT_DIR_EXON
    for SAMPLE in ${SAMPLES}; do
        ${PYTHON_BINARY} ${TOOL_MERGE_FEATURECOUNTS_TABLES} ${SAMPLE}_${PID}_*.featureCounts.dexseq.tsv
    done
fi

if [ "$RUN_LIBRARY_QC" == true ]
then
    mkdir -p ${LIBRARY_QC_DIR}
    for SAMPLE in ${SAMPLES}; do
        ${PYTHON_BINARY} ${TOOL_SINGLE_CELL_LIBRARY_QC} ${SAMPLE} ${PID} ${COUNT_DIR} ${LIBRARY_QC_DIR}
    done
fi