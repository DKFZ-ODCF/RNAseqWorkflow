#!/usr/bin/env bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/RNAseqWorkflow).
#

IFS=$'\n\t'

########################################################################
##                                                                    ##
##                           FUNCTIONS                                ##
##                                                                    ##
########################################################################

 source $TOOL_BASH_LIB

########################################################################
##                                                                    ##
##                           WORK FLOW                                ##
##                                                                    ##
########################################################################


    remove_directory $RODDY_SCRATCH/*
    remove_file $RNASEQC_DIR/${SAMPLE}_${pid}/refGene.txt*
    remove_file $RNASEQC_DIR/${SAMPLE}_${pid}/rRNA_intervals.list
    if [[ -f "$RNASEQC_DIR/${SAMPLE}_${pid}/${SAMPLE}_${pid}_RNAseQC.tgz" ]]
    then
        echo "# RNAseQC results already archived... skipping"
    else
        echo_run "tar --remove-files -czvf $RNASEQC_DIR/${SAMPLE}_${pid}/${SAMPLE}_${pid}_RNAseQC.tgz $RNASEQC_DIR/${SAMPLE}_${pid}/${SAMPLE}_${pid} "
    fi
    remove_file $ALIGNMENT_DIR/$STAR_SORTED_BAM
    remove_file $ALIGNMENT_DIR/$STAR_NOTSORTED_BAM
    remove_file $ALIGNMENT_DIR/$STAR_CHIMERA_SAM
    remove_file $ALIGNMENT_DIR/$STAR_CHIMERA_BAM
	remove_file $ALIGNMENT_DIR/*fifo.read1
	remove_file $ALIGNMENT_DIR/*fifo.read2
	make_directory $ALIGNMENT_DIR/${SAMPLE}_${pid}_star_logs_and_files
	echo_run "mv $ALIGNMENT_DIR/${SAMPLE}_${pid}*out $ALIGNMENT_DIR/${SAMPLE}_${pid}_star_logs_and_files 2>/dev/null"
	echo_run "mv $ALIGNMENT_DIR/${SAMPLE}_${pid}*.tab $ALIGNMENT_DIR/${SAMPLE}_${pid}_star_logs_and_files 2>/dev/null"
	echo_run "mv $ALIGNMENT_DIR/${SAMPLE}_${pid}_merged._STARgenome $ALIGNMENT_DIR/${SAMPLE}_${pid}_star_logs_and_files 2>/dev/null"
	echo_run "mv $ALIGNMENT_DIR/${SAMPLE}_${pid}_merged._STARpass1 $ALIGNMENT_DIR/${SAMPLE}_${pid}_star_logs_and_files 2>/dev/null"
	# remove_directory $SCRATCH


