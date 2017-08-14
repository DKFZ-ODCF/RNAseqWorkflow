/*
 * Copyright (c) 2017 eilslabs.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/Roddy/LICENSE.txt).
 */

package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.files.*
import de.dkfz.roddy.StringConstants
import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.config.ConfigurationValue
import de.dkfz.roddy.core.Analysis
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.core.Workflow
import de.dkfz.roddy.execution.io.BaseMetadataTable
import de.dkfz.roddy.execution.io.MetadataTableFactory
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider
import de.dkfz.roddy.knowledge.files.BaseFile
import de.dkfz.roddy.knowledge.files.FileObject
import de.dkfz.roddy.plugins.LibrariesFactory
import de.dkfz.roddy.tools.Tuple2
import groovy.transform.CompileStatic
import org.apache.commons.io.IOUtils

/**
 * RNA seq workflow based on STAR!
 *
 * Uses some methods from QCPipeline(AlignmentAndQCWorkflows:1.1.51+)
 */
@CompileStatic
class RNAseqWorkflowForSingleCell extends Workflow {

    /**
     * Roddy currently does not have fine grained methods for the MetadataTableFactory.
     * That's why I just extracted it from there. It should be available in Roddy in the future.
     * We mark it deprecated here.
     * @param context
     * @return
     */
    @Deprecated
    Tuple2<Map<String, String>, List<String>> extraceMetadataTableColumnsFromConfig(ExecutionContext context) {
        Analysis analysis = context.analysis
        def missingColValues = []
        def mandatoryColumns = []
        def cvalues = analysis.getConfiguration().getConfigurationValues()

        Map<String, String> columnIDMap = cvalues.get("metadataTableColumnIDs").getValue()
                .split(StringConstants.COMMA)
                .collectEntries {
            String colVar ->
                ConfigurationValue colVal = cvalues.get(colVar);
                if (!colVal) {
                    missingColValues << colVar;
                }

                if (colVal.hasTag("mandatory")) mandatoryColumns << colVal.id;
                return ["${colVar}": colVal?.toString()]
        }
        return new Tuple2(columnIDMap, mandatoryColumns as List<String>);
    }

    @Override
    boolean execute(ExecutionContext context) {
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)

        // Only one sample is allowed, see checks below.
        Sample sample = samples[0]

        RNAseqLaneFileGroupSet lfgs = new RNAseqLaneFileGroupSet(runtimeService.loadLaneFilesForSample(context, sample))

        // Just a dummy, which is needed as a first parameter to some jobs.
        LaneFile dummyFile = lfgs.getFirstLaneFile()

        File barcodeFile = lfgs.getBarcodeFile()

        // Parse welllist using a MetadataTable
        def welllist = getBarcodeFileMDT(context, barcodeFile)

        // Split the table first by its samples, then by the barcodes. This is
        // to get chunks, which are not too large for the Bash command size and for the Star parameters

        // Get the CHUNK_SIZE from the configuration
        final int CHUNK_SIZE = context.getConfiguration().getConfigurationValues().get("chunkSize", "64").toInt()
        // Get the count of chunks, which is the # of record divided by the chunk size
        final int countOfChunks = (int) (welllist.size() / CHUNK_SIZE)

        // Loop over all the chunks, get sub tables by the chunk id
        for (int chunk = 0; chunk < countOfChunks + 1; chunk++) {

            // Create and write chunked barcode file
            createAndWriteChunkedBarcodeFile(context, welllist, chunk)

            // Now comes your workflow. If you got any questions, don't hesitate asking me!
            // Please take a look at Naveeds resimplified RNAseq Workflow

//                String runId = lfgs.getRuns()[0]
//
//                TextFile checkpointfile = (TextFile) call("jemultiplexer", dummyFile, "READS_LEFT=${readsLeft}", "READS_RIGHT=${readsRight}", "BARCODE_FILE=${barcodeFilepath}")
//                checkpointfile = (TextFile) call("trimming", checkpointfile, "BARCODE_FILE=${barcodeFilepath}")
//                checkpointfile = (TextFile) call("singlecellAlignment", checkpointfile, "BARCODE_FILE=${barcodeFilepath}", "RUN_ID=${runId}")
//
//                List<TextFile> checkpointfiles = []
//                TextFile checkpointfile_processing
//                for (int chunk in 1..numJobs) {
//                    checkpointfile_processing = (TextFile) call("rnaseqProcessing", checkpointfile, "CHUNK_INDEX=${chunk}")
//                    checkpointfiles.add(checkpointfile_processing)
//                }
//
//                def checkpointfilegroup = new CheckpointFileGroup(checkpointfiles)
//                call("singlecellPostprocessing", checkpointfile, "BARCODE_FILE=${barcodeFilepath}", checkpointfilegroup)
        }

        return true
    }

    /**
     * As it says, this method will pull a subset from the barcode metadata table and
     * write this to storage.
     * @param barcodeFile
     * @param currentChunk
     * @param context
     * @param chunk
     * @return
     */
    public BaseFile createAndWriteChunkedBarcodeFile(ExecutionContext context, BaseMetadataTable barcodeFile, int currentChunk) {
        COProjectsRuntimeService runtimeService = ((COProjectsRuntimeService) context.runtimeService)
        BaseMetadataTable subsetByChunk = barcodeFile.unsafeSubsetByColumn("Chunk", currentChunk.toString());

        // Assemble the table header for the new barcodeFile files
        final String HEADER = barcodeFile.header.join("\t")

        // Write a new barcodeFile for the next calls.
        BaseFile reducedBarcodeFile = BaseFile.constructSourceFile(
                // Load the class TODO check if this is the right class, replace it with the class needed!
                LibrariesFactory.instance.loadRealOrSyntheticClass("WelllistFile", BaseFile.class.name)
                // Get the new barcodeFile file path
                , new File(runtimeService.getTemporaryDirectory(context), "welllist_${currentChunk}.tsv")
                , context)

        // Get the records for the chunked table
        List<Map<String, String>> records = subsetByChunk.records

        // Get the new file text by:
        // Adding the header line
        // Collect all the record lines by joining the row values with \t
        // Flatten the new file content array and join it by \n
        // Finally write the file content to a remote file inside the temp directory for the workflow run
        String reducedContent = [
                HEADER,
                records.collect {
                    Map<String, String> row -> row.collect { String k, String v -> v }.join("\t")
                }
        ].flatten().join("\n")
        FileSystemAccessProvider.instance.writeTextFile(reducedBarcodeFile.path, reducedContent)
        return reducedBarcodeFile
    }

    /**
     * Load the welllist table, which can or will be chunked later on.
     * @param context
     * @param barcodeFile
     * @return
     */
    BaseMetadataTable getBarcodeFileMDT(ExecutionContext context, File barcodeFile) {
        // We still need to adapt the columns to the actual headers.

        StringBuilder welllistText = new StringBuilder()
        String[] content = FileSystemAccessProvider.instance.loadTextFile(barcodeFile)

        final int CHUNK_SIZE = context.getConfiguration().getConfigurationValues().get("chunkSize", "64").toInt()
        String[] columnHeaders = ("Chunk\t" + content[0]).split("[\t]") // Get the real headers and create the MetadataTable column configuration

        // Prepare the configuration to use the current input for the welllist.
        context.configurationValues.putAll(
                columnHeaders.collectEntries {
                    [it, new ConfigurationValue(context.configuration, it, it, "string", "", ["mandatory"])]
                } as Map<String, ConfigurationValue>
        )
        context.configurationValues.put("metadataTableColumnIDs", columnHeaders.join(","))

        welllistText << "Chunk\t" << content[0] << "\n"// Add the chunk column. This is not in the file.
        for (int i = 1; i < content.size(); i++) { // Add the chunk id and the rest of the row.
            int chunk = (int) (i / CHUNK_SIZE)  // Align to CHUNK_SIZE block index (integer)
            welllistText << "${chunk}\t" << content[i] << "\n"
        }

        String barcodeFileContents = welllistText.toString()
        Tuple2<Map<String, String>, List<String>> mdtSettings = extraceMetadataTableColumnsFromConfig(context)

        return MetadataTableFactory.readTable(new InputStreamReader(IOUtils.toInputStream(barcodeFileContents)), "tsv", mdtSettings.x, mdtSettings.y)
    }

    @Override
    boolean checkExecutability(ExecutionContext context) {
        // TODO Not working for single cell. Should be adapted for this workflow!
        //        boolean result = new RNAseqWorkflow().checkExecutability(context);

        //        if (!result) return false

        // Check if the sample count is right. Must be 1
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)

        if (samples.size() == 0 || !samples || !samples[0]) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There was no sample available for dataset ${context.dataSet.id}"))
            return false
        }

        if (samples.size() > 1) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were too many samples available for dataset ${context.dataSet.id}. Only one is allowed."))
            return false
        }

        // Check welllist and see, if there is more than one entry (in addition to the header)

        RNAseqLaneFileGroupSet lfgs = new RNAseqLaneFileGroupSet(runtimeService.loadLaneFilesForSample(context, samples[0]))

        File barcodeFile = lfgs.getBarcodeFile()

        if (!context.fileIsAccessible(barcodeFile)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("The barcode file ${barcodeFile} is not accessible."))
            return false
        }

        int sizeOfBarcodeFile = FileSystemAccessProvider.instance.loadTextFile(barcodeFile).size()
        if (sizeOfBarcodeFile < 2) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("The barcode file ${barcodeFile} is malformed or empty."))
            return false
        }

        return true
    }
}
