package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.common.ParallelizationHelper
import de.dkfz.b080.co.files.LaneFile
import de.dkfz.b080.co.files.LaneFileGroup
import de.dkfz.b080.co.files.Sample
import de.dkfz.b080.co.files.TextFile
import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.core.Workflow
import de.dkfz.roddy.knowledge.files.FileGroup
import de.dkfz.roddy.knowledge.files.FileObject
import de.dkfz.roddy.knowledge.files.IndexedFileObjects
import groovy.transform.CompileStatic

/**
 * RNA seq workflow based on STAR!
 *
 * Uses some methods from QCPipeline(AlignmentAndQCWorkflows:1.1.51+)
 */
@CompileStatic
class RNAseqWorkflow extends Workflow {

    private void run_rnaseq(RNAseqLaneFileGroupSet lfgs, String sampleName) {
        // The file dummy is used by Roddy. It is not actually needed by the starAlignment job itself
        // but Roddy needs one specific file for every created job.
        LaneFile dummyFile = lfgs.getFirstLaneFile()

        //STAR
        String readsSTARLeft = lfgs.getLeftLaneFilesAsCSVs()
        String readsSTARRight = lfgs.getRightLaneFilesAsCSVs()

        // Kalisto
        String readsKallisto = lfgs.getLaneFilesAlternatingWithSpaceSep()

        // Flowcell ids
        String runIDs = lfgs.getFlowCellIDsWithSpaceSep()

        // Lane ids
        String laneIDs = lfgs.getLaneIDsWithSpaceSep()

        // Read groups per pair with " , " separation ( space comma space )
        String readGroups = lfgs.getBamReadGroupLines()

        call("starAlignment", dummyFile, "SAMPLE=${sampleName}", "READS_STAR_LEFT=${readsSTARLeft}", "READS_STAR_RIGHT=${readsSTARRight}", "READS_KALLISTO=${readsKallisto}", "PARM_RUNIDS=${runIDs}", "PARM_LANEIDS=${laneIDs}", "PARM_READGROUPS=${readGroups}")
    }

    @Override
    boolean execute(ExecutionContext context) {
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()
        boolean runSingleCellDemultiplexing = context.getConfiguration().getConfigurationValues().getBoolean("runSingleCellDemultiplexing", false);
        List<Sample> samples = runtimeService.getSamplesForContext(context)

        for (Sample sample : samples) {
            def laneFilesForSample = runtimeService.loadLaneFilesForSample(context, sample)
            if (!laneFilesForSample)
                continue

            def lfgs = new RNAseqLaneFileGroupSet(laneFilesForSample)

            if (runSingleCellDemultiplexing) {
                // LaneFile barcodeFile = lfgs.getFirstLaneFile().path.getParentFile() //TODO: GET CORRECT BARCODE FILE PATH
                // int numOfCells = barcodeFile.path.readLines().size()

                LaneFile dummyFile = lfgs.getFirstLaneFile()

                List<LaneFileGroup> lfg_list_demultiplexed = []
                for (List<String> l: [lfgs.getLeftLaneFiles(), lfgs.getRightLaneFiles(), lfgs.getLaneIDs(), lfgs.getRuns()].transpose() as List<List<String>>) {
                    FileGroup lanefiles = call("jemultiplexer", dummyFile, "READ_LEFT=${l[0]}", "READ_RIGHT=${l[1]}") as FileGroup
                    lfg_list_demultiplexed.add(new LaneFileGroup(context, l[2], l[3], sample, lanefiles.filesInGroup as List<LaneFile>))
                }
                run_rnaseq(new RNAseqLaneFileGroupSet(lfg_list_demultiplexed), sample.name);
            }
            else run_rnaseq(lfgs, sample.name);
        }

        return true
    }

    @Override
    boolean checkExecutability(ExecutionContext context) {
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()

        /** Check sample availability **/
        List<Sample> samples = runtimeService.getSamplesForContext(context)
        if (samples.size() == 0) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were no samples available."))
            return false
        }

        /** Check fastq availability **/
        boolean hasFiles = false;
        for (Sample sample : samples) {
            def laneFilesForSample = runtimeService.loadLaneFilesForSample(context, sample)
            if (!laneFilesForSample) {
                context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("No fastq files found for sample ${sample.name}"))
                continue
            }
            hasFiles = true
            laneFilesForSample.each {
                LaneFileGroup lfg ->
                    lfg.filesInGroup.each {
                        hasFiles &= true//context.fileIsAccessible(it.path)
                    }
            }
        }

        /** No samples, no files, return false otherwise true **/
        return hasFiles;
    }
}
