package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.files.COFileStageSettings
import de.dkfz.b080.co.files.LaneFile
import de.dkfz.b080.co.files.LaneFileGroup
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.knowledge.files.FileStageSettings
import groovy.transform.CompileStatic

/**
 * Extend the LaneFileGroup list from QCPipeline
 * Created by heinold on 05.12.16.
 */
@CompileStatic
class RNAseqLaneFileGroupSetForPairedEnd {
    protected List<LaneFileGroup> laneFileGroupList = null

    protected Map<LaneFile, LaneFile> laneFiles

    protected ExecutionContext context

    protected RNAseqConfig config

    RNAseqLaneFileGroupSetForPairedEnd(List<LaneFileGroup> laneFileGroupList) {
        this.laneFileGroupList = laneFileGroupList
        this.context = laneFileGroupList[0].executionContext
        this.config = new RNAseqConfig(context)
        laneFiles = laneFileGroupList.collectEntries { LaneFileGroup lfg -> return [lfg.filesInGroup[0], lfg.filesInGroup[1]] }
    }

    LaneFile getFirstLaneFile() {
        return laneFiles.values().first() as LaneFile
    }

    String getLeftLaneFilesAsCSVs() {
        laneFiles.keySet().collect { it.path.absolutePath }.join(",")
    }

    String getRightLaneFilesAsCSVs() {
        laneFiles.values().collect { it.path.absolutePath }.join(",")
    }

    String getLaneFilesAlternatingWithSpaceSep() {
        laneFileGroupList.collect { return "${it.filesInGroup[0].path} ${it.filesInGroup[1].path}" }.join(" ")
    }


    String getFlowCellIDsWithSpaceSep() {
        laneFileGroupList.collect { it.getRun().split("[_]")[-1] }.join(" ")
    }

    String getLaneIDsWithSpaceSep() {
        laneFileGroupList.collect { it.getId() }.join(" ")
    }

    String getFlowCellAndLaneIDsWithSpaceSep() {
        laneFileGroupList.collect { "${it.getRun().split("[_]")[-1]} ${it.getId()}" }.join(" ")
    }

    /**
     * Example for paired end: ID:run150326_D00695_0025_BC6B2MACXX_D2826_GATCAGA_L002 LB:${sample}_${pid}PL:ILLUMINA SM:sample_${sample}_${pid} PU:BC6B2MACXX , ID:run... (space-comma-space separated)
     * @return
     */
    String getBamReadGroupLines() {
        laneFileGroupList.collect { LaneFileGroup lfg ->
            def pid = lfg.filesInGroup.first().getDataSet().getId()
            def sample = lfg.getSample().name
            COFileStageSettings stage = lfg.getFilesInGroup().first().getFileStage() as COFileStageSettings
            [
                    ["ID:" + lfg.getRun(), stage.laneId].join("_"),
                    ["LB:" + sample, pid].join("_"),
                    "PL:ILLUMINA",
                    ["SM:sample", sample, pid].join("_"),
                    "PU:" + lfg.getRun().split("[_]")[-1]
            ].join(" ")
        }.join(" , ")
    }
}
