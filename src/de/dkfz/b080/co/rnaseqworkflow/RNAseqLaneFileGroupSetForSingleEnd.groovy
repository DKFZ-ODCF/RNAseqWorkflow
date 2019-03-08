/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/RNAseqWorkflow).
 */

package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.files.COFileStageSettings
import de.dkfz.b080.co.files.LaneFile
import de.dkfz.b080.co.files.LaneFileGroup
import groovy.transform.CompileStatic

/**
 * Extend the LaneFileGroup list from QCPipeline
 * Created by heinold on 05.12.16.
 */
@CompileStatic
class RNAseqLaneFileGroupSetForSingleEnd extends RNAseqLaneFileGroupSetForPairedEnd {

    protected List<LaneFile> laneFilesList = []

    RNAseqLaneFileGroupSetForSingleEnd(List<LaneFileGroup> laneFileGroupList) {
        super(laneFileGroupList)
        laneFilesList = laneFileGroupList.collect { it.filesInGroup[0] }
    }

    LaneFile getFirstLaneFile() {
        return laneFilesList.first() as LaneFile
    }

    @Override
    String getLeftLaneFilesAsCSVs() {
        laneFilesList.collect { it.path.absolutePath }.join(",")
    }

    @Override
    String getRightLaneFilesAsCSVs() {
        return ""
    }

    String getTrimmedLeftLaneFilesAsCSVs(String trimmingOutputDirectory) {
        laneFileGroupList.collect {
            def it -> it.getRun() + "_" + new File(new File(it.filesInGroup[0].path.parentFile.parentFile, trimmingOutputDirectory), it.filesInGroup[0].path.getName())
        }.join(",")
    }

    String getTrimmedRightLaneFilesAsCSVs(String trimmingOutputDirectory) {
        return ""
    }

    @Override
    String getLaneFilesAlternatingWithSpaceSep() {
        laneFileGroupList.collect { return "${it.filesInGroup[0].path}" }.join(" ")
    }

    @Override
    String getFlowCellIDsWithSpaceSep() {
        laneFileGroupList.collect { it.getRun().split("[_]")[-1] }.join(" ")
    }

    @Override
    String getLaneIDsWithSpaceSep() {
        laneFileGroupList.collect { it.getId() }.join(" ")
    }

    @Override
    String getFlowCellAndLaneIDsWithSpaceSep() {
        laneFileGroupList.collect { "${it.getRun().split("[_]")[-1]} ${it.getId()}" }.join(" ")
    }

    /**
     * Example for single end: ID:run150326_D00695_0025_BC6B2MACXX_D2826_GATCAGA_L002_R1_001 LB:${sample}_${pid}PL:ILLUMINA SM:sample_${sample}_${pid} PU:BC6B2MACXX , ID:run... (space-comma-space separated)
     *
     * This differs from the paired end read group method!
     * lane ids cannot be resolved properly for single end files. So we'll take
     * the whole filename to the suffix and ignore the laneid variable in the file stage.
     * @return
     */
    @Override
    String getBamReadGroupLines() {
        laneFilesList.collect { LaneFile lfg ->
            def pid = lfg.getDataSet().getId()
            def sample = lfg.getSample().name
            def fstage = ((COFileStageSettings)lfg.fileStage)
            COFileStageSettings stage = lfg.getFileStage() as COFileStageSettings
            [
                    ["ID:" + fstage.runID, lfg.path.name.split("[.]")[0]].join("_"),
                    ["LB:" + sample, pid].join("_"),
                    "PL:ILLUMINA",
                    ["SM:sample", sample, pid].join("_"),
                    "PU:" + fstage.runID.split("[_]")[-1]
            ].join(" ")
        }.join(" , ")
    }
}
