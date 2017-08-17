#!/usr/bin/env bash

set -vx

##################################################################
##								##
##  HIPO2 RNAseq workflow (alignment)					##
##  Authors: Naveed Ishaque, Barbara Hutter, Sebastian Uhrig	##
##								##
##################################################################

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
	ARIBA_BINARY=arriba
	ARIBA_READTHROUGH_BINARY=extract_read-through_fusions
	ARIBA_DRAW_FUSIONS=draw_fusions.R
fi

# Check software stack before running workflow
check_executable "$STAR_BINARY"
check_executable "$FEATURECOUNTS_BINARY"
check_executable "$SAMBAMBA_BINARY"
check_executable "$SAMTOOLS_BINARY"
check_executable "$RNASEQC_BINARY"
check_executable "$KALLISTO_BINARY"
check_executable "$QUALIMAP_BINARY"
check_executable "$ARIBA_BINARY"
check_executable "$ARIBA_READTHROUGH_BINARY"
check_executable "$ARIBA_DRAW_FUSIONS"

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

env | sort > $DIR_EXECUTION/${PBS_JOBNAME}.alignment.env_dump.txt

##
## STAR 2-PASS ALIGNMENT: 12 core, 50Gb, 6 hours
##

make_directory $SCRATCH

if [ "$RUN_STAR" == true ]
then
	make_directory $ALIGNMENT_DIR
	cd $ALIGNMENT_DIR
	## the STAR temp directory must not exist, otherwise STAR will fail
	remove_directory ${STAR_TMPDIR}
	if [ "$useSingleEndProcessing" == true ] || [ "$runSingleCellWorkflow" == true ]
	then
	    echo_run "${STAR_BINARY} ${STAR_PARAMS} --readFilesIn ${READS_STAR_LEFT}                     --readFilesCommand ${READ_COMMAND} --outSAMattrRGline ${PARM_READGROUPS}"
	else
	    echo_run "${STAR_BINARY} ${STAR_PARAMS} --readFilesIn ${READS_STAR_LEFT} ${READS_STAR_RIGHT} --readFilesCommand ${READ_COMMAND} --outSAMattrRGline ${PARM_READGROUPS}"
	fi
	check_or_die $STAR_NOTSORTED_BAM alignment
	remove_directory ${STAR_TMPDIR}

    if [ "$runSingleCellWorkflow" == false ]; then
        echo_run "mv ${SAMPLE}_${PID}_merged.Chimeric.out.junction ${SAMPLE}_${PID}_chimeric_merged.junction"
    	check_or_die $STAR_SORTED_BAM alignment
	    check_or_die $STAR_CHIMERA_SAM alignment
        # Currently Below is for ARRIBA - won't be needed for scRNAseq data analysis
        ## BAM-erise and sort chimera file: 1 core, 1 hours, 200mb
        echo_run "$SAMTOOLS_BINARY view -Sbh $STAR_CHIMERA_SAM | $SAMTOOLS_BINARY sort - -o $STAR_CHIMERA_BAM_PREF.bam"
        check_or_die ${STAR_CHIMERA_BAM_PREF}.bam chimera-sam-2-bam
        #echo_run "$SAMBAMBA_BINARY markdup -t 1 -l 0 ${STAR_CHIMERA_BAM_PREF}.bam | $SAMTOOLS_BINARY view -h - | $SAMTOOLS_BINARY view -S -b -@ $CORES > ${STAR_CHIMERA_MKDUP_BAM}"
        echo_run "$SAMBAMBA_BINARY markdup -t $CORES ${STAR_CHIMERA_BAM_PREF}.bam ${STAR_CHIMERA_MKDUP_BAM}"
        check_or_die $STAR_CHIMERA_MKDUP_BAM chimera-post-markdups
        remove_file ${STAR_CHIMERA_BAM_PREF}.bam
        echo_run "$SAMBAMBA_BINARY index -t $CORES $STAR_CHIMERA_MKDUP_BAM"
        check_or_die ${STAR_CHIMERA_MKDUP_BAM}.bai chimera-alignment-index
        echo_run "md5sum $STAR_CHIMERA_MKDUP_BAM | cut -f 1 -d ' ' > $STAR_CHIMERA_MKDUP_BAM.md5"
        check_or_die ${STAR_CHIMERA_MKDUP_BAM}.md5 chimera-alignment-md5sums
    fi
fi

# TODO: Add Salmon for Single Cell Workflow - we won't use Kallisto for scRNAseq

if [ "$RUN_CLEANUP" == true ]
then
    remove_file $ALIGNMENT_DIR/$STAR_CHIMERA_SAM
fi

touch "${CHECKPOINT_ALIGNMENT}"