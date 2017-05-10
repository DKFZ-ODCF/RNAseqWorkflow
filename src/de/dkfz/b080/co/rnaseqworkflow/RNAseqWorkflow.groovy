package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.files.COBaseFile
import de.dkfz.b080.co.files.TextFile
import de.dkfz.b080.co.files.LaneFile
import de.dkfz.b080.co.files.LaneFileGroup
import de.dkfz.b080.co.files.Sample
import de.dkfz.b080.co.files.CheckpointFileGroup
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.core.Workflow
import groovy.transform.CompileStatic

/**
 * RNA seq workflow based on STAR!
 *
 * Uses some methods from QCPipeline(AlignmentAndQCWorkflows:1.1.51+)
 */
@CompileStatic
class RNAseqWorkflow extends Workflow {
    private void run_rnaseq(RNAseqLaneFileGroupSet lfgs, String sampleName, ExecutionContext context) {
        // The file dummy is used by Roddy. It is not actually needed by the starAlignment job itself
        // but Roddy needs one specific file for every created job.

        TextFile dummyFile = (TextFile)COBaseFile.constructManual(TextFile, lfgs.getFirstLaneFile())

        //STAR
        String readsLeft = lfgs.getLeftLaneFilesAsCSVs()
        String readsRight = lfgs.getRightLaneFilesAsCSVs()

        // TODO: BELOW SHOULD BE FIXED
        String readsSTARLeft = lfgs.getTrimmedLeftLaneFilesAsCSVs(context.getConfiguration().getConfigurationValues().getString("trimmingOutputDirectory"))
        String readsSTARRight = lfgs.getTrimmedRightLaneFilesAsCSVs(context.getConfiguration().getConfigurationValues().getString("trimmingOutputDirectory"))

        // Kalisto
        String readsKallisto = lfgs.getLaneFilesAlternatingWithSpaceSep()

        // Flowcell ids
        String runIDs = lfgs.getFlowCellIDsWithSpaceSep()

        // Lane ids
        String laneIDs = lfgs.getLaneIDsWithSpaceSep()

        // Read groups per pair with " , " separation ( space comma space )
        String readGroups = lfgs.getBamReadGroupLines()

        TextFile checkpointfile = (TextFile)call("trimming", dummyFile, "READ1=${readsLeft}", "READ2=${readsRight}")
        checkpointfile = (TextFile)call("starAlignment", checkpointfile, "SAMPLE=${sampleName}", "READS_STAR_LEFT=${readsSTARLeft}", "READS_STAR_RIGHT=${readsSTARRight}", "READS_KALLISTO=${readsKallisto}", "PARM_RUNIDS=${runIDs}", "PARM_LANEIDS=${laneIDs}", "PARM_READGROUPS=${readGroups}")
        call("rnaseqProcessing", checkpointfile, "CHUNK_INDEX=1")
    }

    private void run_singlecell(RNAseqLaneFileGroupSet lfgs, ExecutionContext context) {
        LaneFile dummyFile = lfgs.getFirstLaneFile()

        int numJobs = context.getConfiguration().getConfigurationValues().get("numProcessingJobs").toInt()

        String readsLeft = lfgs.getLeftLaneFilesAsCSVs()
        String readsRight = lfgs.getRightLaneFilesAsCSVs()

        String barcodeFilename = context.getConfiguration().getConfigurationValues().get("rawBarcodeFilename")
        String barcodeFilepath = new File(dummyFile.getPath().getParentFile().getParentFile(), barcodeFilename).absolutePath

        String runId = lfgs.getRuns()[0]

        TextFile checkpointfile = (TextFile)call("jemultiplexer", dummyFile, "READS_LEFT=${readsLeft}", "READS_RIGHT=${readsRight}", "BARCODE_FILE=${barcodeFilepath}")
        checkpointfile = (TextFile)call("trimming", checkpointfile, "BARCODE_FILE=${barcodeFilepath}")
        checkpointfile = (TextFile)call("singlecellAlignment", checkpointfile, "BARCODE_FILE=${barcodeFilepath}", "RUN_ID=${runId}")

        List<TextFile> checkpointfiles = []
        TextFile checkpointfile_processing
        for (int i in 1..numJobs) {
            checkpointfile_processing = (TextFile)call("rnaseqProcessing", checkpointfile, "CHUNK_INDEX=${i}")
            checkpointfiles.add(checkpointfile_processing)
        }

        def checkpointfilegroup = new CheckpointFileGroup(checkpointfiles)
        call("singlecellPostprocessing", checkpointfile, "BARCODE_FILE=${barcodeFilepath}", checkpointfilegroup)
    }

    @Override
    boolean execute(ExecutionContext context) {
        boolean runSingleCellWorkflow = getflag(context, "runSingleCellWorkflow", false)

        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)


        for (Sample sample : samples) {
            def laneFilesForSample = runtimeService.loadLaneFilesForSample(context, sample)
            if (!laneFilesForSample)
                continue

            def lfgs = new RNAseqLaneFileGroupSet(laneFilesForSample)

            if (runSingleCellWorkflow)
                run_singlecell(lfgs, context)
            else
                run_rnaseq(lfgs, sample.name, context)
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
        boolean hasFiles = false
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
                        hasFiles &= context.fileIsAccessible(it.path)
                    }
            }
        }

        if (getflag(context, "runSingleCellDemultiplexing", false)) {
            //TODO Introduce the check at the right position
//            context.fileIsAccessible(getBarcodeFile(context, ))
        }

        /** No samples, no files, return false otherwise true **/
        return hasFiles
    }
}
