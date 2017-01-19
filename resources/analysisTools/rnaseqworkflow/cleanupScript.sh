#!/usr/bin/env bash

IFS=$'\n\t'

########################################################################
##                                                                    ##
##                           FUNCTIONS                                ##
##                                                                    ##
########################################################################

 source $TOOL_NAV_LIB

########################################################################
##                                                                    ##
##                           WORK FLOW                                ##
##                                                                    ##
########################################################################


    remove_directory $SCRATCH/*
    remove_file $RNASEQC_DIR/${SAMPLE}_${PID}/refGene.txt*
    remove_file $RNASEQC_DIR/${SAMPLE}_${PID}/rRNA_intervals.list
    if [[ -f "$RNASEQC_DIR/${SAMPLE}_${PID}/${SAMPLE}_${PID}_RNAseQC.tgz" ]]
    then
        echo "# RNAseQC results already archived... skipping"
    else
        echo_run "tar --remove-files -czvf $RNASEQC_DIR/${SAMPLE}_${PID}/${SAMPLE}_${PID}_RNAseQC.tgz $RNASEQC_DIR/${SAMPLE}_${PID}/${SAMPLE}_${PID} "
    fi
    remove_file $ALIGNMENT_DIR/$STAR_SORTED_BAM
    remove_file $ALIGNMENT_DIR/$STAR_NOTSORTED_BAM
    remove_file $ALIGNMENT_DIR/$STAR_CHIMERA_SAM
    remove_file $ALIGNMENT_DIR/$STAR_CHIMERA_BAM
	remove_file $ALIGNMENT_DIR/*fifo.read1
	remove_file $ALIGNMENT_DIR/*fifo.read2
	make_directory $ALIGNMENT_DIR/${SAMPLE}_${PID}_star_logs_and_files
	echo_run "mv $ALIGNMENT_DIR/${SAMPLE}_${PID}*out $ALIGNMENT_DIR/${SAMPLE}_${PID}_star_logs_and_files 2>/dev/null"
	echo_run "mv $ALIGNMENT_DIR/${SAMPLE}_${PID}*.tab $ALIGNMENT_DIR/${SAMPLE}_${PID}_star_logs_and_files 2>/dev/null"
	echo_run "mv $ALIGNMENT_DIR/${SAMPLE}_${PID}_merged._STARgenome $ALIGNMENT_DIR/${SAMPLE}_${PID}_star_logs_and_files 2>/dev/null"
	echo_run "mv $ALIGNMENT_DIR/${SAMPLE}_${PID}_merged._STARpass1 $ALIGNMENT_DIR/${SAMPLE}_${PID}_star_logs_and_files 2>/dev/null"
	remove_directory $SCRATCH


