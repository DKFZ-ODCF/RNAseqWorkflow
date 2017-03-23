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
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider
import de.dkfz.roddy.knowledge.files.FileGroup
import de.dkfz.roddy.knowledge.files.FileObject
import de.dkfz.roddy.knowledge.files.IndexedFileObjects
import groovy.transform.CompileStatic
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser

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

    private boolean run_singlecell(RNAseqLaneFileGroupSet lfgs, ExecutionContext context, Sample sample) {
        LaneFile dummyFile = lfgs.getFirstLaneFile()
        File barcodeFile = getBarcodeFile(context, dummyFile)
        def barcodeFileContent = FileSystemAccessProvider.getInstance().loadTextFile(barcodeFile)

        // Read stuff from csv, tsv   convertFormat is in Roddy! Take Roddy from GitHub
        //  => look up BaseMetadataTableFactory
        //CSVFormat tableFormat = convertFormat(format)
        //tableFormat = tableFormat.withCommentMarker('#' as char)
        //        .withIgnoreEmptyLines()
        //        .withHeader();
        //CSVParser parser = tableFormat.parse(instream)
        int numOfCells = barcodeFileContent.size()

        if (numOfCells == 0) return false;

        List<LaneFileGroup> lfg_list_demultiplexed = []
        List<List<String>> ll = [lfgs.getLeftLaneFiles(), lfgs.getRightLaneFiles(), lfgs.getLaneIDs(), lfgs.getRuns()].transpose()

        for (List<String> l in ll) {
            // It is not written wrong: It reall is called jemultiplefacetxer! From JE demultiplexer
            FileGroup lanefiles = callWithOutputFileGroup("jemultiplexer", dummyFile, numOfCells, "READ_LEFT=${l[0]}", "READ_RIGHT=${l[1]}") as FileGroup
            lfg_list_demultiplexed.add(new LaneFileGroup(context, l[2], l[3], sample, lanefiles.filesInGroup as List<LaneFile>))
        }
        run_rnaseq(new RNAseqLaneFileGroupSet(lfg_list_demultiplexed), sample.name)
    }

    @Override
    boolean execute(ExecutionContext context) {
        boolean runSingleCellDemultiplexing = getflag(context, "runSingleCellDemultiplexing", false)

        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)

        for (Sample sample : samples) {
            def laneFilesForSample = runtimeService.loadLaneFilesForSample(context, sample)
            if (!laneFilesForSample)
                continue

            def lfgs = new RNAseqLaneFileGroupSet(laneFilesForSample)

            if (runSingleCellDemultiplexing) {
                return run_singlecell(lfgs, context, sample)
            } else run_rnaseq(lfgs, sample.name)
        }

        return true
    }

    private File getBarcodeFile(ExecutionContext context, LaneFile dummyFile) {
        String barcodeFilename = context.getConfiguration().getConfigurationValues().get("barcodeFilename")
        File barcodeFile = new File(dummyFile.getPath().getParentFile(), barcodeFilename)
        return barcodeFile
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
