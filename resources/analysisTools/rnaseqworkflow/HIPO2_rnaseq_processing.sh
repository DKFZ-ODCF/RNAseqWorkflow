#!/usr/bin/env bash

set -vx

# TODO: I don't like below two lines
SAMPLES=( `cat ${CHECKPOINT_ALIGNMENT} | cut -f 1` )
STAR_SORTED_BAMS=( `cat ${CHECKPOINT_ALIGNMENT} | cut -f 2` )

total_cnt="${#STAR_SORTED_BAMS[@]}"
let chunk_len="(${total_cnt}+${numProcessingJobs}-1)/${numProcessingJobs}"
let start="(${CHUNK_INDEX}-1)*${chunk_len}"

SAMPLES=( ${SAMPLES[@]:${start}:${chunk_len}} )
STAR_SORTED_BAMS=( ${STAR_SORTED_BAMS[@]:${start}:${chunk_len}} )

if [[ ${#STAR_SORTED_BAMS[@]} -eq 0  ]]; then
    echo "There is no BAM file to process."
    touch ${CHECKPOINT_PROCESSING}
    exit 0
fi

if [ "$runSingleCellWorkflow" == true ]; then
   	STAR_SORTED_MKDUP_BAMS=()
    for STAR_SORTED_BAM in "${STAR_SORTED_BAMS[@]}"; do
        STAR_SORTED_MKDUP_BAM=${STAR_SORTED_BAM%.bam}.mdup.bam
        STAR_SORTED_MKDUP_BAMS+=(${STAR_SORTED_MKDUP_BAM})
    done
else
    STAR_SORTED_MKDUP_BAMS=${STAR_SORTED_MKDUP_BAM}
fi

##################################################################
##								##
##  HIPO2 RNAseq workflow (processing)					##
##  Authors: Naveed Ishaque, Barbara Hutter, Sebastian Uhrig	##
##								##
##################################################################

# TODO : 21/03/2017 : a : (a.1) parse read length, (a.2) specificy sbjdOverhang based on read length and (a.2) pick appropiate genome index
# TODO : 21/03/2017 : b : (b.1) parse strand specificity from RNAseQC, (b.2) only compute read counts for the correct direction using feature counts

##################################################################
##								##
##			   SETUP ENV				##
##								##
##################################################################

source $CONFIG_FILE

source $TOOL_NAV_LIB

echo_run $DO_FIRST

##
## SOFTWARE STACK : load module # module load HIPO2_rna/v1 # loads default software versions for the HIPO2 RNAseq workflow
##

if [ "$LOAD_MODULE" == true ]
then
	module load $MODULE_ENV

	STAR_BINARY=STAR
	FEATURECOUNTS_BINARY=featureCounts
	SAMBAMBA_BINARY=sambamba
	SAMTOOLS_BINARY=samtools # samtools_bin defaults to v0.0.19 despite the RNAseq.XML!!!!!!!!!!!!!!!!!!!!!!!!!!
	RNASEQC_BINARY=rnaseqc
	KALLISTO_BINARY=kallisto
	QUALIMAP_BINARY=qualimap
	ARRIBA_BINARY=arriba
	ARRIBA_READTHROUGH_BINARY=extract_read-through_fusions
	ARRIBA_DRAW_FUSIONS=draw_fusions.R
fi

# Check software stack before running workflow
check_executable "$STAR_BINARY"
check_executable "$FEATURECOUNTS_BINARY"
check_executable "$SAMBAMBA_BINARY"
check_executable "$SAMTOOLS_BINARY"
check_executable "$RNASEQC_BINARY"
check_executable "$KALLISTO_BINARY"
check_executable "$QUALIMAP_BINARY"
check_executable "$ARRIBA_BINARY"
check_executable "$ARRIBA_READTHROUGH_BINARY"
check_executable "$ARRIBA_DRAW_FUSIONS"

########################################################################
##								    ##
##			   WORK FLOW				##
##								    ##
########################################################################

set -u

#make_directory $RESULTS_PID_DIR

##
## PRINT env
##

env | sort > $DIR_EXECUTION/${PBS_JOBNAME}.processing.env_dump.txt

make_directory $SCRATCH

if [ "$RUN_STAR" == true ]
then
    make_directory $ALIGNMENT_DIR
    cd $ALIGNMENT_DIR

    ## markdups using sambamba  (requires 7Gb and 20 min walltime (or 1.5 hrs CPU time) for 200m reads)
    parallel_run 2 STAR_SORTED_BAMS STAR_SORTED_MKDUP_BAMS ${CORES} "${SAMBAMBA_BINARY} markdup -t1 \$0 \$1"
    for STAR_SORTED_MKDUP_BAM in "${STAR_SORTED_MKDUP_BAMS[@]}"; do
        check_or_die $STAR_SORTED_MKDUP_BAM post-markdups
    done

    ## index using samtools (requires 40MB and 5 minutes for 200m reads)
    parallel_run 1 STAR_SORTED_MKDUP_BAMS ${CORES} "${SAMBAMBA_BINARY} index -t1 \$0"
    for STAR_SORTED_MKDUP_BAM in "${STAR_SORTED_MKDUP_BAMS[@]}"; do
        check_or_die $STAR_SORTED_MKDUP_BAM.bai alignment-index
    done

    ## md5sum
    parallel_run 1 STAR_SORTED_MKDUP_BAMS ${CORES} "md5sum \$0 | cut -f 1 -d \" \"  > \$0.md5"
    for STAR_SORTED_MKDUP_BAM in "${STAR_SORTED_MKDUP_BAMS[@]}"; do
        check_or_die $STAR_SORTED_MKDUP_BAM.md5 alignment-md5sums
    done

    ## flagstats (requires 4MB and 5 minutes for 200m reads)
    parallel_run 1 STAR_SORTED_MKDUP_BAMS ${CORES} "$SAMBAMBA_BINARY flagstat -t 1 \$0 > \$0.flagstat"
    for STAR_SORTED_MKDUP_BAM in "${STAR_SORTED_MKDUP_BAMS[@]}"; do
        check_or_die $STAR_SORTED_MKDUP_BAM.flagstat alignment-qc-flagstats
    done
fi

# Run the fingerprinting. This requires the .bai file.
if [[ "${runFingerprinting:-false}" == true ]]
then
    cd $ALIGNMENT_DIR
    parallel_run 1 STAR_SORTED_MKDUP_BAMS ${CORES} "$TOOL_FINGERPRINT $fingerprintingSitesFile \$0 > \$0.fp.tmp && mv \$0.fp.tmp \$0.fp"
fi

##
## RNAseQC (requires 10Gb and 6 hours for 200m reads)
##

if [ "$RUN_RNASEQC" == true ]
then
    make_directory $RNASEQC_DIR/${PID}
    cd $RNASEQC_DIR/${PID}
    DOC_FLAG=" "
    if [ "$disableDoC_GATK" == true ]
    then
		DOC_FLAG="-noDoC"
	fi

    echo_run paste <(printf "%s\n" ${SAMPLES[@]}) <(printf "%s\n" ${STAR_SORTED_MKDUP_BAMS[@]}) | awk -v PID=$PID '{split($2, s, "_"); print $1 "_" PID "\t" $1 "\t" s[4] "_" s[5] "_" s[6]}' > rnaseqc_samples.txt
    echo_run "$RNASEQC_BINARY -r $GENOME_GATK_INDEX $DOC_FLAG -t $GENE_MODELS -n 1000 -o . -s rnaseqc_samples.txt &> $DIR_EXECUTION/${PBS_JOBNAME}.${PID}_RNAseQC.log &"
fi

##
## QualiMap2 (requires 6Gb and 3 hours for 200m reads)
##

if [ "$RUN_QUALIMAP" == true ]
then
	make_directory $QUALIMAP_DIR/${SAMPLES}_${PID} # TODO: THIS IS NOT VALID FOR SINGLE CELL WORKFLOW
	cd $QUALIMAP_DIR/${SAMPLES}_${PID}
    for STAR_SORTED_MKDUP_BAM in "${STAR_SORTED_MKDUP_BAMS[@]}"; do
    	if [ "$runSingleCellWorkflow" == true ]; then
	        QUALIMAP_OUTFILE=${STAR_SORTED_MKDUP_BAM}.report
	        QUALIMAP_BAM=$ALIGNMENT_DIR/$STAR_SORTED_MKDUP_BAM
    	else
	        QUALIMAP_OUTFILE=${SAMPLES}_${PID}.report
	        QUALIMAP_BAM=$ALIGNMENT_DIR/$STAR_NOTSORTED_BAM
    	fi
        echo_run "$QUALIMAP_BINARY rnaseq -gtf $GENE_MODELS -s -pe --java-mem-size=60G -outfile ${QUALIMAP_OUTFILE} -outdir $QUALIMAP_DIR/${SAMPLES}_${PID} -bam ${QUALIMAP_BAM}"
   	done
    check_or_die rnaseq_qc_results.txt qc-qualimap2
fi

##
## feature Counts htseq like requires 1gb, and 30 min of 4 cores
##

if [ "$RUN_FEATURE_COUNTS" == true ]
then
	make_directory $COUNT_DIR 
	cd $COUNT_DIR
	COUNT="-t exon -g gene_id -Q 255 -T $CORES -a $GENE_MODELS -F GTF --donotsort"

    if [ "$runSingleCellWorkflow" == true ]; then
        FEATURECOUNTS_PREFIXES=(${STAR_SORTED_MKDUP_BAMS[@]})
        FEATURECOUNTS_BAMS=(${STAR_SORTED_MKDUP_BAMS[@]})
    else
        FEATURECOUNTS_PREFIXES=${SAMPLES}_${PID}
        FEATURECOUNTS_BAMS=$STAR_NOTSORTED_BAM
    	COUNT="${COUNT} -p -B"
    fi

    for S in {0..2}
    do
        parallel_run 2 FEATURECOUNTS_PREFIXES FEATURECOUNTS_BAMS ${CORES} "mkdir -p $SCRATCH/featureCounts_\$0_$S && $FEATURECOUNTS_BINARY $COUNT --tmpDir $SCRATCH/featureCounts_\$0_$S -s $S -o \$0.featureCounts.s$S $ALIGNMENT_DIR/\$1 && rm -fr $SCRATCH/featureCounts_\$0_$S"
        for FEATURECOUNTS_PREFIX in "${FEATURECOUNTS_PREFIXES[@]}"; do
            check_or_die ${FEATURECOUNTS_PREFIX}.featureCounts.s${S} gene-counting
        done
    done

    # RPKM TPM calculations & cleanup
    parallel_run 1 FEATURECOUNTS_PREFIXES ${CORES} "$TOOL_COUNTS_TO_FPKM_TPM \$0.featureCounts.s0 \$0.featureCounts.s1 \$0.featureCounts.s2 $GENE_MODELS $GENE_MODELS_EXCLUDE > \$0.fpkm_tpm.featureCounts.tsv"
    if [ "$runSingleCellWorkflow" == true ]; then
        rm -f non_mapped_reads.txt
    fi
    for FEATURECOUNTS_PREFIX in "${FEATURECOUNTS_PREFIXES[@]}"; do
    	check_or_die ${FEATURECOUNTS_PREFIX}.fpkm_tpm.featureCounts.tsv counting-featureCounts
        if [ "$runSingleCellWorkflow" == true ]; then
        	echo_run "echo -e \"${FEATURECOUNTS_PREFIX}\t`grep NoFeatures ${FEATURECOUNTS_PREFIX}.featureCounts.s0.summary | cut -f 2`\" >> non_mapped_reads.txt"
        fi
    	make_directory ${FEATURECOUNTS_PREFIX}_featureCounts_raw
        echo_run "mv ${FEATURECOUNTS_PREFIX}.featureCounts* ${FEATURECOUNTS_PREFIX}_featureCounts_raw"
    done

	parallel_run 1 FEATURECOUNTS_PREFIXES ${CORES} "tar --remove-files -czvf \$0_featureCounts_raw.tgz \$0_featureCounts_raw"
    #if [ "$runSingleCellWorkflow" == true ]; then
    #    echo_run ${PYTHON_BINARY} ${TOOL_MERGE_FEATURECOUNTS_TABLES} ${SAMPLES}_${PID}_*.featureCounts.tsv
    #fi
fi

##
## feature Counts DEXSEQ like
##

if [ "$RUN_FEATURE_COUNTS_DEXSEQ" == true ]
then
	make_directory $COUNT_DIR_EXON
	cd $COUNT_DIR_EXON
	COUNT_EXONS="-f -O -F GTF -a $GENE_MODELS_DEXSEQ -t exonic_part -p -Q 255 -T $CORES --tmpDir $SCRATCH/${SAMPLE}_${pid}_featureCountsExons --donotsort"
	make_directory $SCRATCH/${SAMPLES}_${PID}_featureCountsExons

    if [ "$runSingleCellWorkflow" == true ]; then
        FEATURECOUNTS_DEXSEQ_PREFIXES=(${STAR_SORTED_MKDUP_BAMS[@]})
        FEATURECOUNTS_DEXSEQ_BAMS=(${STAR_SORTED_MKDUP_BAMS[@]})
    else
        FEATURECOUNTS_DEXSEQ_PREFIXES=${SAMPLES}_${PID}
        FEATURECOUNTS_DEXSEQ_BAMS=$STAR_NOTSORTED_BAM
    fi

    for S in {0..2}
    do
        parallel_run 2 FEATURECOUNTS_DEXSEQ_PREFIXES FEATURECOUNTS_DEXSEQ_BAMS ${CORES} "mkdir -p $SCRATCH/featureCounts_\$0_$S && $FEATURECOUNTS_BINARY $COUNT_EXONS --tmpDir $SCRATCH/featureCounts_\$0_$S -s $S -o \$0.featureCounts.s$S $ALIGNMENT_DIR/\$1 && rm -fr $SCRATCH/featureCounts_\$0_$S"
        for FEATURECOUNTS_DEXSEQ_PREFIX in "${FEATURECOUNTS_DEXSEQ_PREFIXES[@]}"; do
    		check_or_die ${FEATURECOUNTS_DEXSEQ_PREFIX}.featureCounts.dexseq.s${S} exon-counting
    	done
    done

    ## RPKM TPM calculations
    parallel_run 1 FEATURECOUNTS_DEXSEQ_PREFIXES ${CORES} "$TOOL_COUNTSDEXSEQ_TO_FPKM_TPM \$0.featureCounts.s0 \$0.featureCounts.s1 \$0.featureCounts.s2 $GENE_MODELS $GENE_MODELS_EXCLUDE > \$0.fpkm_tpm.featureCounts.dexseq.tsv"

    for FEATURECOUNTS_DEXSEQ_PREFIX in "${FEATURECOUNTS_DEXSEQ_PREFIXES[@]}"; do
	    check_or_die ${FEATURECOUNTS_DEXSEQ_PREFIX}.fpkm_tpm.featureCounts.dexseq.tsv counting-featureCounts_dexseq
    	# cleanup
	    make_directory ${FEATURECOUNTS_DEXSEQ_PREFIX}_featureCounts_dexseq_raw
    	echo_run "mv ${FEATURECOUNTS_DEXSEQ_PREFIX}.featureCounts* ${FEATURECOUNTS_DEXSEQ_PREFIX}_featureCounts_dexseq_raw"
	done
    parallel_run 1 FEATURECOUNTS_DEXSEQ_PREFIXES ${CORES} "tar --remove-files -czvf \$0_featureCounts_dexseq_raw.tgz \$0_featureCounts_dexseq_raw"
	#remove_directory $SCRATCH/${SAMPLE}_${pid}_featureCountsExons
fi

# TODO: Array conversion is needed to store large number of Fastq files
# https://stackoverflow.com/questions/14525296/bash-check-if-variable-is-array
#if [[ !("$(declare -p READS_KALLISTO)" =~ "declare -a") ]]; then
#    READS_KALLISTO=( $READS_KALLISTO ) # Make sure it is array
#fi

# TODO: BELOW KALLISTO COMMANDS ARE NOT SO SUITABLE FOR SINGLE CELL WORKFLOW
KALLISTO_PARAMS="quant -i $GENOME_KALLISTO_INDEX -o . -t $CORES -b 100"
if [ "$RUN_KALLISTO" == true ]
then
	make_directory $KALLISTO_UN_DIR/${SAMPLES}_${pid}
	cd $KALLISTO_UN_DIR/${SAMPLES}_${pid}
	echo_run "$KALLISTO_BINARY $KALLISTO_PARAMS $READS_KALLISTO"
	check_or_die abundance.tsv kallisto
	echo_run "$TOOL_KALLISTO_RESCALE abundance.tsv $GENE_MODELS_EXCLUDE > abundance.rescaled.tsv"
fi

if [ "$RUN_KALLISTO_RF" == true ]
then
	make_directory $KALLISTO_RF_DIR
	cd $KALLISTO_RF_DIR
	$KALLISTO_BINARY $KALLISTO_PARAMS --rf-stranded $READS_KALLISTO
	check_or_die abundance.tsv kallisto
	$TOOL_KALLISTO_RESCALE abundance.tsv $GENE_MODELS_EXCLUDE > abundance.rescaled.tsv
fi

if [ "$RUN_KALLISTO_FR" == true ]
then
	make_directory $KALLISTO_FR_DIR
	cd $KALLISTO_FR_DIR
	$KALLISTO_BINARY $KALLISTO_PARAMS --fr-stranded $READS_KALLISTO
	check_or_die abundance.tsv kallisto
	$TOOL_KALLISTO_RESCALE abundance.tsv $GENE_MODELS_EXCLUDE > abundance.rescaled.tsv
fi

##
## Fusion detection
##

if [ "$RUN_ARRIBA" == true ]
then
	make_directory $ARRIBA_DIR
	cd $ARRIBA_DIR

    for STAR_SORTED_MKDUP_BAM in "${STAR_SORTED_MKDUP_BAMS[@]}"; do
        if [ "$runSingleCellWorkflow" == true ]; then
            ARRIBA_PREFIX=${STAR_SORTED_MKDUP_BAM%_merged.mdup.bam}
        else
            ARRIBA_PREFIX=${SAMPLES}_${PID}
        fi
	    echo_run "$ARRIBA_READTHROUGH_BINARY -g $GENE_MODELS -i $ALIGNMENT_DIR/$STAR_SORTED_MKDUP_BAM -o ${ARRIBA_PREFIX}_merged_read_through.bam"
	    echo_run "$ARRIBA_BINARY -c $ALIGNMENT_DIR/$STAR_CHIMERA_MKDUP_BAM -r ${ARRIBA_PREFIX}_merged_read_through.bam -x $ALIGNMENT_DIR/$STAR_SORTED_MKDUP_BAM -a $GENOME_FA -k $ARRIBA_KNOWN_FUSIONS -g $GENE_MODELS -b $ARRIBA_BLACKLIST -o ${ARRIBA_PREFIX}.fusions.txt -O ${ARRIBA_PREFIX}.discarded_fusions.txt "
	    if [[ -f "${ARRIBA_PREFIX}.fusions.txt" ]]
	    then
		    echo_run "$ARRIBA_DRAW_FUSIONS --annotation=$GENE_MODELS --fusions=${ARRIBA_PREFIX}.fusions.txt --output=${ARRIBA_PREFIX}.fusions.pdf"
	    fi
	done
fi

##
## Produce json QC file
##

if [ "$RUN_QCJSON" == true ]
then
	make_directory $QC_DIR
	if [ "$RUN_RNASEQC" == true ]
	then
		echo "# Waiting for RNAseQC"
		wait && cd $QC_DIR
	else
		cd $QC_DIR
	fi

	if [[ -f "$RNASEQC_DIR/${SAMPLES}_${pid}/metrics.tsv" ]]
	then
		echo_run "mv $RNASEQC_DIR/${SAMPLES}_${pid}/metrics.tsv $RNASEQC_DIR/${SAMPLES}_${pid}/${SAMPLES}_${pid}_metrics.tsv"
	fi

	if  [[ -f "$RNASEQC_DIR/${SAMPLES}_${pid}/${SAMPLES}_${pid}_metrics.tsv" ]]
	then
		echo "#FOUND FILE: RNAseQC file \"$RNASEQC_DIR/${SAMPLES}_${pid}/${SAMPLES}_${pid}_metrics.tsv\""
	else
		echo "#ERROR: file not found: \"$RNASEQC_DIR/${SAMPLES}_${pid}/${SAMPLES}_${pid}_metrics.tsv\" ... exitting!" 1>&2
		exit 1
	fi

	if [[ -f "$DIR_EXECUTION/${PBS_JOBNAME}.${SAMPLES}_${pid}_RNAseQC.log" ]]
	then
		check_text_and_die $DIR_EXECUTION/${PBS_JOBNAME}.${SAMPLES}_${pid}_RNAseQC.log "org.broadinstitute.sting.gatk.walkers.coverage.DepthOfCoverageWalker.onTraversalDone" "rerun with \$disableDoC_GATK=TRUE"
		check_text_or_die $DIR_EXECUTION/${PBS_JOBNAME}.${SAMPLES}_${pid}_RNAseQC.log "Finished Successfully"
	fi
	echo_run "$TOOL_CREATE_JSON_FROM_OUTPUT $ALIGNMENT_DIR/${STAR_SORTED_MKDUP_BAM}.flagstat $RNASEQC_DIR/${SAMPLES}_${pid}/${SAMPLES}_${pid}_metrics.tsv > ${JSON_PREFIX}qualitycontrol.json"
	check_or_die ${JSON_PREFIX}qualitycontrol.json qc-json

# TODO: ABOVE SHOULD BE MODIFIED LIKE BELOW
#    for STAR_SORTED_MKDUP_BAM in "${STAR_SORTED_MKDUP_BAMS[@]}"; do
#        if [ "$runSingleCellWorkflow" == true ]; then
#            QCJSON_PREFIX=${STAR_SORTED_MKDUP_BAM%_merged.mdup.bam}
#            RNASEQC_DIR_PER_SAMPLE=$RNASEQC_DIR/${QCJSON_PREFIX}
#            JSON_PREFIX=${JSON_PREFIX}${QCJSON_PREFIX}
#        else
#            QCJSON_PREFIX=${SAMPLES}_${PID}
#            RNASEQC_DIR_PER_SAMPLE=$RNASEQC_DIR/${QCJSON_PREFIX}
#        fi
#
#    	if [[ -f "${RNASEQC_DIR_PER_SAMPLE}/metrics.tsv" ]]
#	    then
#    		echo_run "mv ${RNASEQC_DIR_PER_SAMPLE}/metrics.tsv ${RNASEQC_DIR_PER_SAMPLE}/${QCJSON_PREFIX}_metrics.tsv"
#	    fi
#	    echo_run "$TOOL_CREATE_JSON_FROM_OUTPUT $ALIGNMENT_DIR/${STAR_SORTED_MKDUP_BAM}.flagstat ${RNASEQC_DIR_PER_SAMPLE}/${QCJSON_PREFIX}_metrics.tsv > ${JSON_PREFIX}qualitycontrol.json"
#	    check_or_die ${JSON_PREFIX}qualitycontrol.json qc-json
#	done
fi

##
## Clean up files which are no longer needed
##

if [ "$RUN_CLEANUP" == true ]
then
	if [ "$RUN_RNASEQC" == true ]
	then
		echo "# Waiting for RNAseQC"
		wait
	fi
	remove_directory $SCRATCH/*

    for STAR_SORTED_MKDUP_BAM in "${STAR_SORTED_MKDUP_BAMS[@]}"; do
        if [ "$runSingleCellWorkflow" == true ]; then
            RNASEQC_DIR_PER_SAMPLE=$RNASEQC_DIR/${STAR_SORTED_MKDUP_BAM%_merged.mdup.bam}
        else
            RNASEQC_DIR_PER_SAMPLE=$RNASEQC_DIR/${PID}
        fi

        if [[ -f "${RNASEQC_DIR_PER_SAMPLE}/${PID}_RNAseQC.tgz" ]]
	    then
	        echo "# RNAseQC results already archived... skipping"
    	else
	        echo_run "tar --remove-files -czvf ${RNASEQC_DIR_PER_SAMPLE}/${PID}_RNAseQC.tgz ${RNASEQC_DIR_PER_SAMPLE}/${PID} "
	    fi

    	remove_file ${RNASEQC_DIR_PER_SAMPLE}/refGene.txt*
	    remove_file ${RNASEQC_DIR_PER_SAMPLE}/rRNA_intervals.list
	done

	for STAR_SORTED_BAM in "${STAR_SORTED_BAMS[@]}"; do
	    remove_file $ALIGNMENT_DIR/$STAR_SORTED_BAM
	done
	remove_file $ALIGNMENT_DIR/$STAR_NOTSORTED_BAM
	remove_file $ALIGNMENT_DIR/$STAR_CHIMERA_SAM
	remove_file $ALIGNMENT_DIR/*fifo.read1
	remove_file $ALIGNMENT_DIR/*fifo.read2
	make_directory $ALIGNMENT_DIR/${PID}_star_logs_and_files
	echo_run "mv -f $ALIGNMENT_DIR/${PID}*out $ALIGNMENT_DIR/${PID}_star_logs_and_files 2>/dev/null"
	echo_run "mv -f $ALIGNMENT_DIR/${PID}*.tab $ALIGNMENT_DIR/${PID}_star_logs_and_files 2>/dev/null"
	echo_run "mv -f $ALIGNMENT_DIR/${PID}_merged._STARgenome $ALIGNMENT_DIR/${PID}_star_logs_and_files 2>/dev/null"
	echo_run "mv -f $ALIGNMENT_DIR/${PID}_merged._STARpass1 $ALIGNMENT_DIR/${PID}_star_logs_and_files 2>/dev/null"
	remove_directory $SCRATCH
fi

touch ${CHECKPOINT_PROCESSING}

echo "DONE!"
