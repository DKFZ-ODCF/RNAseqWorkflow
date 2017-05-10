export SAMPLE

source $CONFIG_FILE

SAMPLES=`${PYTHON_BINARY} ${TOOL_PARSE_BARCODE_FILE} -t2 ${BARCODE_FILE}`

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
        ${RSCRIPT_BINARY} ${TOOL_SINGLE_CELL_LIBRARY_QC} ${COUNT_DIR}/${SAMPLE}_${PID}_featureCounts.count.tsv ${TOOL_MAKE_LIBRARY_QC_PLOT} ${LIBRARY_QC_DIR} ${SAMPLE}
    done
fi