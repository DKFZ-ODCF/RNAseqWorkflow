package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.files.COFileStageSettings
import de.dkfz.b080.co.files.LaneFile
import de.dkfz.b080.co.files.LaneFileGroup
import de.dkfz.b080.co.files.Sample
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
        laneFiles = laneFileGroupList.collectEntries { LaneFileGroup lfg -> return [lfg.filesInGroup[0], (lfg.filesInGroup.size() >= 2) ? lfg.filesInGroup[1] : null] }
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
        laneFileGroupList.collect {
            def it -> it.getRun() + "_" + new File(new File(it.filesInGroup[0].path.parentFile.parentFile, trimmingOutputDirectory), it.filesInGroup[0].path.getName())
        }.join(",")
    }

    String getTrimmedRightLaneFilesAsCSVs(String trimmingOutputDirectory) {
        laneFileGroupList.collect {
            def it -> it.getRun() + "_" + new File(new File(it.filesInGroup[1].path.parentFile.parentFile, trimmingOutputDirectory), it.filesInGroup[1].path.getName())
        }.join(",")
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

    String getPID() {
        laneFileGroupList.first().filesInGroup.first().getDataSet().getId()
    }

    /**
     * Get the barcode file (aka welllist) from the configuration.
     * @return
     */
    File getBarcodeFile() {
        ExecutionContext context = getFirstLaneFile().executionContext
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()
        Sample sample = ((COFileStageSettings)getFirstLaneFile().fileStage).sample

        RNAseqLaneFileGroupSetForPairedEnd lfgs = new RNAseqLaneFileGroupSetForPairedEnd(runtimeService.loadLaneFilesForSample(context, sample))
        LaneFile dummyFile = lfgs.getFirstLaneFile()

        String barcodeFilename = context.getConfiguration().getConfigurationValues().get("rawBarcodeFilename")
        File barcodeFile = new File(dummyFile.getPath().getParentFile().getParentFile(), barcodeFilename)
        return barcodeFile
    }

    Map<String, List<LaneFile>> getSingleCellLaneFilesByRun() {
        def runs = laneFileGroupList.collect { it.run }.sort().unique()
        Map<String, List<LaneFile>> result = [:]

        runs.each { String run ->
            result[run] =
                    laneFileGroupList.findAll { LaneFileGroup lfg ->
                        lfg.run == run
                    }.collect {
                        LaneFileGroup group ->
                            group.filesInGroup[0] // For single cell, the first one is always the valid file. The second one is a dummy
                    } as List<LaneFile>
        }
        return result
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
