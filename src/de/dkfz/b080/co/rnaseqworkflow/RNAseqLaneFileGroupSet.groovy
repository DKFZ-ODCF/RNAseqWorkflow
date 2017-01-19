package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.files.COFileStageSettings
import de.dkfz.b080.co.files.LaneFile
import de.dkfz.b080.co.files.LaneFileGroup
import de.dkfz.roddy.knowledge.files.FileStageSettings
import groovy.transform.CompileStatic

/**
 * Extend the LaneFileGroup list from QCPipeline
 * Created by heinold on 05.12.16.
 */
@CompileStatic
class RNAseqLaneFileGroupSet {
    private List<LaneFileGroup> laneFileGroupList = null;

    private Map<LaneFile, LaneFile> laneFiles

    RNAseqLaneFileGroupSet(List<LaneFileGroup> laneFileGroupList) {
        this.laneFileGroupList = laneFileGroupList
        laneFiles = laneFileGroupList.collectEntries { LaneFileGroup lfg -> return [lfg.filesInGroup[0], lfg.filesInGroup[1]] }
    }

    public LaneFile getFirstLaneFile() {
        return laneFiles?.values()?.first() as LaneFile
    }

    public String getLeftLaneFilesAsCSVs() {
        laneFiles.keySet().collect { it.path.absolutePath }.join(",");
    }

    public String getRightLaneFilesAsCSVs() {
        laneFiles.values().collect { it.path.absolutePath }.join(",");
    }

    public String getLaneFilesAlternatingWithSpaceSep() {
        laneFileGroupList.collect { return "${it.filesInGroup[0].path} ${it.filesInGroup[1].path}" }.join(" ")
    }


    public String collectFlowCellIDs() {
        laneFileGroupList.collect { it.getRun().split("[_]")[-1] }.join(" ")
    }

    public String collectLaneIDs() {
        laneFileGroupList.collect { it.getId() }.join(" ")
    }

    public String collectFlowCellAndLaneIDsWithSpaceSep() {
        laneFileGroupList.collect { "${it.getRun().split("[_]")[-1]} ${it.getId()}" }.join(" ")
    }

    public String collectReadGroups() {
        //    ID:run150326_D00695_0025_BC6B2MACXX_${PID}_D2826_GATCAGA_L002 LB:${sample}_${pid}PL:ILLUMINA SM:sample_${sample}_${pid} PU:BC6B2MACXX , ID:run... (space-comma-space separated)

        laneFileGroupList.collect { LaneFileGroup lfg ->
            def pid = lfg.filesInGroup.first().getDataSet().getId()
            def sample = lfg.getSample().name
            COFileStageSettings stage = lfg.getFilesInGroup().first().getFileStage() as COFileStageSettings
            [
                    ["ID:" + lfg.getRun(), pid, stage.laneId].join("_"),
                    ["LB:" + sample, pid].join("_"),
                    "PL:ILLUMINA",
                    ["SM:" + sample, pid].join("_"),
                    "PU:" + lfg.getRun().split("[_]")[-1]
            ].join(" ")
        }.join(" , ")
    }

}
