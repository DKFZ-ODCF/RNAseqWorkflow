package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.files.COBaseFile
import de.dkfz.b080.co.files.TextFile
import de.dkfz.b080.co.files.LaneFile
import de.dkfz.b080.co.files.LaneFileGroup
import de.dkfz.b080.co.files.Sample
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

    @Override
    boolean execute(ExecutionContext context) {
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()
        RNAseqConfig config = new RNAseqConfig(context)

        List<Sample> samples = runtimeService.getSamplesForContext(context)

        // The file dummy is used by Roddy. It is not actually needed by the starAlignment job itself
        // but Roddy needs one specific file for every created job.
        for (Sample sample : samples) {
            def laneFilesForSample = runtimeService.loadLaneFilesForSample(context, sample)
            if (!laneFilesForSample)
                continue

            def lfgs

            if (config.useSingleEndProcessing()) {
                lfgs = new RNAseqLaneFileGroupSetForSingleEnd(laneFilesForSample)
            } else {
                lfgs = new RNAseqLaneFileGroupSetForPairedEnd(laneFilesForSample)
            }

            // The file dummy is used by Roddy. It is not actually needed by the starAlignment job itself
            // but Roddy needs one specific file for every created job to create dependencies between jobs.
            TextFile dummyFile = (TextFile) COBaseFile.constructManual(TextFile, lfgs.getFirstLaneFile())

            //STAR
            String readsLeft = lfgs.getLeftLaneFilesAsCSVs()
            String readsRight = lfgs.getRightLaneFilesAsCSVs()

            // TODO: BELOW SHOULD BE FIXED
            String readsSTARLeft = lfgs.getTrimmedLeftLaneFilesAsCSVs(context.getConfiguration().getConfigurationValues().getString("trimmingOutputDirectory"))
            String readsSTARRight = config.usePairedEndProcessing() ? lfgs.getTrimmedRightLaneFilesAsCSVs(context.getConfiguration().getConfigurationValues().getString("trimmingOutputDirectory")) : "n/a"

            // Kalisto
            String readsKallisto = lfgs.getLaneFilesAlternatingWithSpaceSep()

            // Read groups per pair with " , " separation ( space comma space )
            String readGroups = lfgs.getBamReadGroupLines()

            TextFile checkpointfile = (TextFile) call("trimming", dummyFile, "SAMPLE=${sample.name}", "READ1=${readsLeft}", "READ2=${readsRight}")
            call("rnaseqProcessing", checkpointfile, "SAMPLE=${sample.name}", "READS_STAR_LEFT=${readsSTARLeft}", "READS_STAR_RIGHT=${readsSTARRight}", "READS_KALLISTO=${readsKallisto}", "PARM_READGROUPS=${readGroups}")
        }
        return true;
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
                    if (new RNAseqConfig(context).usePairedEndProcessing()) {
                        lfg.filesInGroup.each {
                            hasFiles &= context.fileIsAccessible(it.path)
                        }
                    } else { // Single end. Just check the first file.
                        hasFiles &= context.fileIsAccessible(lfg.filesInGroup[0].path)
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
