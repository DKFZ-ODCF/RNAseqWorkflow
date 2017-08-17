/*
 * Copyright (c) 2017 eilslabs.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/Roddy/LICENSE.txt).
 */

package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.common.MetadataTable
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
        String pid = lfgs.getLaneIDs().first().split("[_]").first()
        String runid = lfgs.getRuns().first()

        // Parse welllist using a MetadataTable
        Map<String, Map<String, Map<String, BaseMetadataTable>>> barcodeMDTBySampleAndChunk = getChunkedBarcodeFileMDTablesBySample(context, barcodeFile)

        // Split the table first by its samples, then by the barcodes. This is
        // to get chunks, which are not too large for the Bash command size and for the Star parameters

        String barcodeFilePath = barcodeFile.absolutePath

        String readsLeft = lfgs.getLeftLaneFilesAsCSVs()
        String readsRight = lfgs.getRightLaneFilesAsCSVs()

        BaseFile jeInputFile = createAndWriteJeInputFile(context, pid, barcodeMDTBySampleAndChunk)
        String jeInputFilePath = jeInputFile.absolutePath
        //TextFile checkpointfile_jemultiplexer = (TextFile) call("jemultiplexer", dummyFile, "READS_LEFT=${readsLeft}", "READS_RIGHT=${readsRight}", "JE_INPUT_FILE=${jeInputFilePath}")

        TextFile checkpointfile;
        for (String sampleID : barcodeMDTBySampleAndChunk.keySet()) {
            List<TextFile> checkpointfiles_alignment = []
            for(String chunkID : barcodeMDTBySampleAndChunk[sampleID]["alignment"].keySet()) {
                BaseFile alnInputFile = createAndWriteAlignmentInputFile(context, pid, runid, chunkID, barcodeMDTBySampleAndChunk[sampleID]["alignment"][chunkID])
                //barcodeFilePath = barcodeFileBySampleAndChunk.absolutePath
                //checkpointfile = (TextFile) call("trimming", checkpointfile_jemultiplexer, "BARCODE_FILE=${barcodeFilePath}", "CHUNK_INDEX=${chunkID}")
                //checkpointfile = (TextFile) call("singlecellAlignment", checkpointfile, "BARCODE_FILE=${barcodeFilePath}", "SAMPLE=${sampleID}", "CHUNK_INDEX=${chunkID}")
                //checkpointfiles_alignment.add(checkpointfile)
            }
            //def checkpointfilegroup_alignment = new CheckpointFileGroup(checkpointfiles_alignment)
            List<TextFile> checkpointfiles_processing = []
            for(String chunkID : barcodeMDTBySampleAndChunk[sampleID]["processing"].keySet()) {
                BaseFile barcodeFileBySampleAndChunk = createAndWriteChunkedBarcodeFile(context, barcodeMDTBySampleAndChunk[sampleID]["alignment"][chunkID], sampleID, "alignment", chunkID)
                //barcodeFilePath = barcodeFileBySampleAndChunk.absolutePath
                //checkpointfile = (TextFile) call("rnaseqProcessing", dummyFile, "BARCODE_FILE=${barcodeFilePath}", "CHUNK_INDEX=${chunkID}", checkpointfilegroup_alignment)
                //checkpointfiles_processing.add(checkpointfile)
            }
            def checkpointfilegroup_processing = new CheckpointFileGroup(checkpointfiles_processing)
            BaseFile barcodeFileBySample = createAndWriteChunkedBarcodeFile(context, barcodeMDTBySampleAndChunk[sampleID]["postprocessing"]["all"], sampleID, "postprocessing", "all")
            //barcodeFilePath = barcodeFileBySample.absolutePath
            call("singlecellPostprocessing", dummyFile, "BARCODE_FILE=${barcodeFilePath}", checkpointfilegroup_processing)
        }
        return true
    }

    public BaseFile createAndWriteAlignmentInputFile(ExecutionContext context, String pid, String runid, String chunkID, BaseMetadataTable barcodeTable) {
        // Creates input fastq file, read group, BAM file names
        COProjectsRuntimeService runtimeService = ((COProjectsRuntimeService) context.runtimeService)
        List<Map<String, String>> records = barcodeTable.records
        BaseFile inputFile = BaseFile.constructSourceFile(
                LibrariesFactory.instance.loadRealOrSyntheticClass("TextFile", BaseFile.class.name)
                , new File(runtimeService.getTemporaryDirectory(context), "alignment_"+chunkID+".in")
                , context)
        String fileContent = [
                records.collect {
                    Map<String, String> row -> [
                            pid+"-"+row["Sample"]+"-C"+row["_SampleNum"]+"-H"+sprintf('%03d', row["_CellNum"].toInteger())+"_"+row["Barcode"]+"_R2.fastq.gz",
                            "ID:"+runid+"_C"+row["_SampleNum"]+"_H"+sprintf('%03d', row["_CellNum"].toInteger())+"_"+row["Barcode"]+" LB:"+row["Sample"]+"_"+pid+" PL:ILLUMINA SM:sample_"+row["Sample"]+"_"+pid+" PU:"+runid.split("[_]")[-1]
                            row["Sample"]+"_"+pid+"_C"+row["_SampleNum"]+"_H"+sprintf('%03d', row["_CellNum"].toInteger())+"_"+row["Barcode"]+".demultiplexed.bam"
                    ].join("\t")
                }
        ].flatten().join("\n") + "\n"
        FileSystemAccessProvider.instance.writeTextFile(inputFile.path, fileContent)
        return inputFile
    }

    public BaseFile createAndWriteJeInputFile(ExecutionContext context, String pid, Map<String, Map<String, Map<String, BaseMetadataTable>>> barcodeTables) {
        COProjectsRuntimeService runtimeService = ((COProjectsRuntimeService) context.runtimeService)
        BaseFile inputFile = BaseFile.constructSourceFile(
                LibrariesFactory.instance.loadRealOrSyntheticClass("TextFile", BaseFile.class.name)
                , new File(runtimeService.getTemporaryDirectory(context), "jemultiplexer.in")
                , context)
        String fileContent = ""
        for (String sampleID : barcodeTables.keySet()) {
            BaseMetadataTable barcodeTable = barcodeTables[sampleID]["postprocessing"]["all"]
            List<Map<String, String>> records = barcodeTable.records
            fileContent += [
                    records.collect {
                        Map<String, String> row -> [
                                row["Sample"]+"-R"+row["Row"]+"C"+row["Col"],
                                row["Barcode"],
                                pid+"-"+row["Sample"]+"-C"+row["_SampleNum"]+"-H"+sprintf('%03d', row["_CellNum"].toInteger())+"_"+row["Barcode"]+"_R1.fastq.gz",
                                pid+"-"+row["Sample"]+"-C"+row["_SampleNum"]+"-H"+sprintf('%03d', row["_CellNum"].toInteger())+"_"+row["Barcode"]+"_R2.fastq.gz"
                        ].join("\t")
                    }
            ].flatten().join("\n") + "\n"
        }
        FileSystemAccessProvider.instance.writeTextFile(inputFile.path, fileContent)
        return inputFile
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
    public BaseFile createAndWriteChunkedBarcodeFile(ExecutionContext context, BaseMetadataTable barcodeFile, String sampleID, String typename, String chunkID) {
        COProjectsRuntimeService runtimeService = ((COProjectsRuntimeService) context.runtimeService)

        // Assemble the table header for the new barcodeFile files
        final String HEADER = barcodeFile.header.join("\t")

        // Write a new barcodeFile for the next calls.
        BaseFile reducedBarcodeFile = BaseFile.constructSourceFile(
                // Load the class TODO check if this is the right class, replace it with the class needed!
                LibrariesFactory.instance.loadRealOrSyntheticClass("WelllistFile", BaseFile.class.name)
                // Get the new barcodeFile file path
                , new File(runtimeService.getTemporaryDirectory(context), "welllist_${sampleID}_${typename}_${chunkID}.tsv")
                , context)

        // Get the records for the chunked table
        List<Map<String, String>> records = barcodeFile.records

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
    Map<String, Map<String, Map<String, BaseMetadataTable>>> getChunkedBarcodeFileMDTablesBySample(ExecutionContext context, File barcodeFile) {
        // We still need to adapt the columns to the actual headers.

        StringBuilder welllistText = new StringBuilder()
        String[] content = FileSystemAccessProvider.instance.loadTextFile(barcodeFile)

        final int CHUNK_SIZE = context.getConfiguration().getConfigurationValues().get("chunkSize", "512").toInt()
        final int PROCESSING_JOBS = 20
        String[] columnHeaders = ("_AlignmentChunk\t_ProcessingChunk\t_SampleNum\t_CellNum\t" + content[0]).split("[\t]") // Get the real headers and create the MetadataTable column configuration

        // Prepare the configuration to use the current input for the welllist.
        context.configurationValues.putAll(
                columnHeaders.collectEntries {
                    [it, new ConfigurationValue(context.configuration, it, it, "string", "", ["mandatory"])]
                } as Map<String, ConfigurationValue>
        )
        context.configurationValues.put("metadataTableColumnIDs", columnHeaders.join(","))

        welllistText << "_AlignmentChunk\t_ProcessingChunk\t_SampleNum\t_CellNum\t" << content[0] << "\n"// Add the chunk column. This is not in the file.
        for (int i = 1; i < content.size(); i++) { // Add zero as the chunk id
            welllistText << "0\t0\t0\t0\t" << content[i] << "\n"
        }

        String barcodeFileContents = welllistText.toString()
        Tuple2<Map<String, String>, List<String>> mdtSettings = extraceMetadataTableColumnsFromConfig(context)

        BaseMetadataTable fullTable = MetadataTableFactory.readTable(new InputStreamReader(IOUtils.toInputStream(barcodeFileContents)), "tsv", mdtSettings.x, mdtSettings.y)

        //Modify the records
        Map<String, Map<String, Map<String, BaseMetadataTable>>> allTablesBySampleAndChunk = [:]
        // Get all samples, then all sample tables, then set their chunk ids.
        List<String> sampleIdentifiers = fullTable.listColumn("Sample").sort().unique()
        int sampleNum = 1
        for (String sample : sampleIdentifiers) {
            int cellNum = 1
            BaseMetadataTable tableBySample = fullTable.unsafeSubsetByColumn("Sample", sample)
            allTablesBySampleAndChunk[sample] = [:]
            // Chunk for STAR / Salmon
            for (int i = 0; i < tableBySample.size(); i++) {
                int chunkID = (int)(i / CHUNK_SIZE)
                tableBySample.records[i]["_AlignmentChunk"] = chunkID.toString()
                tableBySample.records[i]["_SampleNum"] = sampleNum.toString()
                tableBySample.records[i]["_CellNum"] = (cellNum++).toString()
            }
            // Chunk for Processing
            int numChunks = PROCESSING_JOBS;
            if (tableBySample.size() < PROCESSING_JOBS) {
                numChunks = tableBySample.size()
            } else if (((int)(Math.ceil(tableBySample.size() / (float)PROCESSING_JOBS))) > CHUNK_SIZE) {
                numChunks = (int)(tableBySample.size() / CHUNK_SIZE)
            }
            int chunkSize = (int)(Math.ceil(tableBySample.size() / (float)numChunks))
            for (int i = 0; i < tableBySample.size(); i++) {
                String chunkID = ((int)(i / chunkSize) + 1).toString()
                tableBySample.records[i]["_Chunk"] = chunkID
            }
            // Store metadata tables
            allTablesBySampleAndChunk[sample]["postprocessing"] = [:]
            allTablesBySampleAndChunk[sample]["postprocessing"]["all"] = tableBySample
            allTablesBySampleAndChunk[sample]["alignment"] = [:]
            for (String chunk in tableBySample.listColumn("_AlignmentChunk").sort().unique() ) {
                allTablesBySampleAndChunk[sample]["alignment"][chunk] = tableBySample.unsafeSubsetByColumn("_AlignmentChunk", chunk)
            }
            allTablesBySampleAndChunk[sample]["processing"] = [:]
            for (String chunk in tableBySample.listColumn("_ProcessingChunk").sort().unique() ) {
                allTablesBySampleAndChunk[sample]["processing"][chunk] = tableBySample.unsafeSubsetByColumn("_ProcessingChunk", chunk)
            }
            sampleNum++
        }
        return allTablesBySampleAndChunk
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
