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
        String pid = lfgs.getPID()
        String runid = lfgs.getRuns().first()

        // Parse welllist using a MetadataTable
        Map<String, Map<String, BaseMetadataTable>> barcodeMDTBySampleAndChunk = getChunkedBarcodeFileMDTablesBySample(context, barcodeFile)

        // Split the table first by its samples, then by the barcodes. This is
        // to get chunks, which are not too large for the Bash command size and for the Star parameters

        String barcodeFilePath = barcodeFile.absolutePath

        String readsLeft = lfgs.getLeftLaneFilesAsCSVs()
        String readsRight = lfgs.getRightLaneFilesAsCSVs()

        BaseFile jeInputFile = createAndWriteJeInputFile(context, pid, barcodeMDTBySampleAndChunk)
        String jeInputFilePath = jeInputFile.absolutePath
        TextFile checkpointfile_jemultiplexer = (TextFile) call("jemultiplexer", dummyFile, "READS_LEFT=${readsLeft}", "READS_RIGHT=${readsRight}", "BARCODE_JE=${jeInputFilePath}")

        TextFile checkpointfile;
        for (String sampleID : barcodeMDTBySampleAndChunk.keySet()) {
            List<TextFile> checkpointfiles_alignment = []
            for(String chunkID : barcodeMDTBySampleAndChunk[sampleID].keySet()) {
                if (chunkID.equals("all"))
                    continue
                List<String> strs = createInputFileStringAndReadGroups(context, pid, runid, chunkID, barcodeMDTBySampleAndChunk[sampleID][chunkID])
                checkpointfile = (TextFile) call("trimming", checkpointfile_jemultiplexer, "READ1=${strs[0]}", "SAMPLE=${sampleID}", "CHUNK_INDEX=${chunkID}")
                checkpointfile = (TextFile) call("rnaseqProcessing", checkpointfile, "SAMPLE=${sampleID}", "READS_STAR_LEFT=${strs[1]}", "PARM_READGROUPS=${strs[2]}","CHUNK_INDEX=${chunkID}")
                checkpointfiles_alignment.add(checkpointfile)
            }
            int numChunks = barcodeMDTBySampleAndChunk[sampleID].keySet().size() - 1
            def checkpointfilegroup_alignment = new CheckpointFileGroup(checkpointfiles_alignment)
            call("singlecellPostprocessing", dummyFile, "BAM_FILE_NUM=${numChunks}", "SAMPLE=${sampleID}", checkpointfilegroup_alignment)
        }
        return true
    }

    public List<String> createInputFileStringAndReadGroups(ExecutionContext context, String pid, String runid, String chunkID, BaseMetadataTable barcodeTable) {
        // Creates input fastq file, read group
        COProjectsRuntimeService runtimeService = ((COProjectsRuntimeService) context.runtimeService)
        String outdir = context.getOutputDirectory().absolutePath
        String demultiplexdir = context.getConfiguration().getConfigurationValues().getString("demultiplexOutputDirectory")
        String trimmingdir = context.getConfiguration().getConfigurationValues().getString("trimmingOutputDirectory")
        List<Map<String, String>> records = barcodeTable.records
        List<List<String>> rtn = [
                records.collect {
                    Map<String, String> row ->
                        [
                            outdir + "/" + demultiplexdir + "/" + pid + "-" + row["Sample"] + "-C" + row["_SampleNum"] + "-H" + sprintf('%03d', row["_CellNum"].toInteger()) + "_" + row["Barcode"] + "_R2.fastq.gz",
                            outdir + "/" + trimmingdir + "/" + pid + "-" + row["Sample"] + "-C" + row["_SampleNum"] + "-H" + sprintf('%03d', row["_CellNum"].toInteger()) + "_" + row["Barcode"] + "_R2.fastq.gz",
                            "ID:" + runid + "_C" + row["_SampleNum"] + "_H" + sprintf('%03d', row["_CellNum"].toInteger()) + "_" + row["Barcode"] + " LB:" + row["Sample"] + "_" + pid + " PL:ILLUMINA SM:sample_" + row["Sample"] + "_" + pid + " PU:" + (runid.split("[_]")[-1])
                        ]
                }
        ][0].transpose()
        return [rtn[0].join(" "), rtn[1].join(","), rtn[2].join(" , ")]
    }

    public BaseFile createAndWriteJeInputFile(ExecutionContext context, String pid, Map<String, Map<String, BaseMetadataTable>> barcodeTables) {
        COProjectsRuntimeService runtimeService = ((COProjectsRuntimeService) context.runtimeService)
        BaseFile inputFile = BaseFile.constructSourceFile(
                LibrariesFactory.instance.loadRealOrSyntheticClass("TextFile", BaseFile.class.name)
                , new File(runtimeService.getTemporaryDirectory(context), "jemultiplexer.in")
                , context)
        String fileContent = ""
        for (String sampleID : barcodeTables.keySet()) {
            BaseMetadataTable barcodeTable = barcodeTables[sampleID]["all"]
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
     * Load the welllist table, which can or will be chunked later on.
     * @param context
     * @param barcodeFile
     * @return
     */
    Map<String, Map<String, BaseMetadataTable>> getChunkedBarcodeFileMDTablesBySample(ExecutionContext context, File barcodeFile) {
        // We still need to adapt the columns to the actual headers.

        StringBuilder welllistText = new StringBuilder()
        String[] content = FileSystemAccessProvider.instance.loadTextFile(barcodeFile)

        final int CHUNK_SIZE = context.getConfiguration().getConfigurationValues().get("chunkSize", "512").toInt()
        String[] columnHeaders = ("_Chunk\t_SampleNum\t_CellNum\t" + content[0]).split("[\t]") // Get the real headers and create the MetadataTable column configuration

        // Prepare the configuration to use the current input for the welllist.
        context.configurationValues.putAll(
                columnHeaders.collectEntries {
                    [it, new ConfigurationValue(context.configuration, it, it, "string", "", ["mandatory"])]
                } as Map<String, ConfigurationValue>
        )
        context.configurationValues.put("metadataTableColumnIDs", columnHeaders.join(","))

        welllistText << "_Chunk\t_SampleNum\t_CellNum\t" << content[0] << "\n"// Add the chunk column. This is not in the file.
        for (int i = 1; i < content.size(); i++) { // Add zero as the chunk id
            welllistText << "0\t0\t0\t" << content[i] << "\n"
        }

        String barcodeFileContents = welllistText.toString()
        Tuple2<Map<String, String>, List<String>> mdtSettings = extraceMetadataTableColumnsFromConfig(context)

        BaseMetadataTable fullTable = MetadataTableFactory.readTable(new InputStreamReader(IOUtils.toInputStream(barcodeFileContents)), "tsv", mdtSettings.x, mdtSettings.y)

        //Modify the records
        Map<String, Map<String, BaseMetadataTable>> allTablesBySampleAndChunk = [:]
        // Get all samples, then all sample tables, then set their chunk ids.
        List<String> sampleIdentifiers = fullTable.listColumn("Sample").sort().unique()
        int sampleNum = 1
        for (String sample : sampleIdentifiers) {
            int cellNum = 1
            BaseMetadataTable tableBySample = fullTable.unsafeSubsetByColumn("Sample", sample)
            sample = sample.replaceAll(/\s+/, "").replaceAll(/[_-]+/, "")
            allTablesBySampleAndChunk[sample] = [:]
            // Chunk for STAR / Salmon
            for (int i = 0; i < tableBySample.size(); i++) {
                int chunkID = (int)(i / CHUNK_SIZE)
                tableBySample.records[i]["_Chunk"] = chunkID.toString()
                tableBySample.records[i]["_SampleNum"] = sampleNum.toString()
                tableBySample.records[i]["_CellNum"] = (cellNum++).toString()
                tableBySample.records[i]["Sample"] = sample
            }
            // Store metadata tables
            allTablesBySampleAndChunk[sample]["all"] = tableBySample
            for (String chunk in tableBySample.listColumn("_Chunk").sort().unique()) {
                allTablesBySampleAndChunk[sample][chunk] = tableBySample.unsafeSubsetByColumn("_Chunk", chunk)
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
