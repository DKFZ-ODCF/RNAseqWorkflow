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
class RNAseqLaneFileGroupSet {
    private List<LaneFileGroup> laneFileGroupList = null

    private Map<LaneFile, LaneFile> laneFiles

    RNAseqLaneFileGroupSet(List<LaneFileGroup> laneFileGroupList) {
        this.laneFileGroupList = laneFileGroupList
        laneFiles = laneFileGroupList.collectEntries { LaneFileGroup lfg -> return [lfg.filesInGroup[0], (lfg.filesInGroup.size()>=2)?lfg.filesInGroup[1]:null] }
    }

    List<String> getLeftLaneFiles() {
        laneFiles.keySet().collect { it.path.absolutePath }
    }

    List<String> getRightLaneFiles() {
        laneFiles.values().collect { it.path.absolutePath }
    }

    List<String> getRuns() {
        laneFileGroupList.collect { it.getRun() }
    }

    List<String> getLaneIDs() {
        laneFileGroupList.collect { it.getId() }
    }

    LaneFile getFirstLaneFile() {
        return laneFiles.keySet().first() as LaneFile
    }

    String getLeftLaneFilesAsCSVs() {
        getLeftLaneFiles().join(",")
    }

    String getRightLaneFilesAsCSVs() {
        getRightLaneFiles().join(",")
    }

    String getTrimmedLeftLaneFilesAsCSVs(String trimmingOutputDirectory) {
        laneFiles.keySet().collect { new File(new File(it.path.parentFile.parentFile, trimmingOutputDirectory), it.path.getName()) }.join(",")
    }

    String getTrimmedRightLaneFilesAsCSVs(String trimmingOutputDirectory) {
        laneFiles.values().collect { new File(new File(it.path.parentFile.parentFile, trimmingOutputDirectory), it.path.getName()) }.join(",")
    }

    String getLaneFilesAlternatingWithSpaceSep() {
        laneFileGroupList.collect { return "${it.filesInGroup[0].path} ${it.filesInGroup[1].path}" }.join(" ")
    }

    String getFlowCellIDsWithSpaceSep() {
        laneFileGroupList.collect { it.getRun().split("[_]")[-1] }.join(" ")
    }

    String getLaneIDsWithSpaceSep() {
        getLaneIDs().join(" ")
    }

    String getFlowCellAndLaneIDsWithSpaceSep() {
        laneFileGroupList.collect { "${it.getRun().split("[_]")[-1]} ${it.getId()}" }.join(" ")
    }

    /**
     * Example: ID:run150326_D00695_0025_BC6B2MACXX_D2826_GATCAGA_L002 LB:${sample}_${pid} PL:ILLUMINA SM:sample_${sample}_${pid} PU:BC6B2MACXX , ID:run... (space-comma-space separated)
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
