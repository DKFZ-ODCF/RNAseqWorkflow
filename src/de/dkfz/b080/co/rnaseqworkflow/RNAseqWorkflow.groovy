package de.dkfz.b080.co.rnaseqworkflow;

import de.dkfz.b080.co.common.BasicCOProjectsRuntimeService;
import de.dkfz.b080.co.common.COProjectsRuntimeService;
import de.dkfz.b080.co.common.WorkflowUsingMergedBams;
import de.dkfz.b080.co.files.BasicBamFile;
import de.dkfz.b080.co.files.COConstants
import de.dkfz.b080.co.files.LaneFile
import de.dkfz.b080.co.files.LaneFileGroup;
import de.dkfz.b080.co.files.Sample;
import de.dkfz.b080.co.qcworkflow.QCPipeline;
import de.dkfz.roddy.config.Configuration;
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues;
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError;
import de.dkfz.roddy.core.Workflow;
import de.dkfz.roddy.knowledge.methods.GenericMethod
import groovy.transform.CompileStatic

import java.lang.reflect.Method;
import java.util.List;

/**
 * RNA seq workflow based on STAR!
 *
 * Uses some methods from QCPipeline(AlignmentAndQCWorkflows:1.1.51+)
 */
@CompileStatic
public class RNAseqWorkflow extends Workflow {

    @Override
    public boolean execute(ExecutionContext context) {
        Configuration cfg = context.getConfiguration();
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService();

        List<Sample> samples = runtimeService.getSamplesForContext(context);
        if (samples.size() == 0)
            return false;

        for (Sample sample : samples) {
            def lfgs = new RNAseqLaneFileGroupSet(runtimeService.loadLaneFilesForSample(context, sample))
            if (!lfgs) {
                context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("No fastq files found for sample ${sample.name}"));
                continue
            };

            LaneFile fakeFile = lfgs.getFirstLaneFile()

            //STAR
            String readsSTARLeft = lfgs.getLeftLaneFilesAsCSVs()
            String readsSTARRight = lfgs.getRightLaneFilesAsCSVs()

            // Kalisto
            String readsKallisto = lfgs.getLaneFilesAlternatingWithSpaceSep()

            // Flowcell ids
            String runIDs = lfgs.collectFlowCellIDs()

            // Lane ids
            String laneIDs = lfgs.collectLaneIDs()

            // Read groups per pair with " , " separation ( space comma space )
            String readGroups = lfgs.collectReadGroups()

            call("starAlignment", fakeFile, "SAMPLE=${sample.name}" , "READS_STAR_LEFT=${readsSTARLeft}", "READS_STAR_RIGHT=${readsSTARRight}", "READS_KALLISTO=${readsKallisto}", "PARM_RUNIDS=${runIDs}", "PARM_LANEIDS=${laneIDs}", "PARM_READGROUPS=${readGroups}")
        }

        return true;
    }
}
