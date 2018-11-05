/*
 * Copyright (c) 2017 eilslabs.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/Roddy/LICENSE.txt).
 */

package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.files.COConstants
import de.dkfz.b080.co.files.Sample
import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.config.ConfigurationConstants
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ContextResource
import de.dkfz.roddy.plugins.LibrariesFactory
import groovy.transform.CompileStatic
import org.junit.BeforeClass
import org.junit.ClassRule
import org.junit.Ignore
import org.junit.Test

import static de.dkfz.b080.co.rnaseqworkflow.Helper.setPrivateField

/**
 * Created by heinold on 05.12.16.
 */
@CompileStatic
class RNAseqLaneFileGroupSetTestForSingleEnd {

    @ClassRule
    public static final ContextResource contextResource = new ContextResource()

    static ExecutionContext context
    static Sample sample0 = new Sample(context, "tumor0")
    static Sample sample1 = new Sample(context, "tumor02")
    static Map<Sample, RNAseqLaneFileGroupSetForPairedEnd> fileGroupSetMap = [:]

    // The test directory
    static String userDir = new File(System.getProperty("user.dir")).parent
    static String singleEndFolderSample0 = userDir + "/out/test/RNAseqWorkflow/resources/testdata_singleEnd/TEST_PID/tumor0/paired/"
    String f0l = singleEndFolderSample0 + "run160319_D00133_0107_BC5YE7ACXX/sequence/D2826_GATCAGA_L002_R1_001.fastq.gz"
    String f1l = singleEndFolderSample0 + "run160326_D00695_0025_BC6B2MACXX/sequence/D2826_GATCAGA_L002_R1_001.fastq.gz"

    @BeforeClass
    static void setup() {
        context = contextResource.createSimpleContext(RNAseqLaneFileGroupSetTestForSingleEnd, new Configuration(null), new COProjectsRuntimeService())

        // The setup is a bit complicated, because I want to load the fastq files from disk to simulate the full
        // fastq loading behaviour

        def resourceDirectory = LibrariesFactory.getGroovyClassLoader().getResource("resources/testdata_singleEnd").file

        setPrivateField("inputDirectory", context, new File(resourceDirectory))

        RNAseqLaneFileGroupSetTestForPairedEnd.setupServices()

        def values = context.getConfiguration().getConfigurationValues()
        values.put("useSingleEndProcessing", "true")
        values.put(ConfigurationConstants.CFG_INPUT_BASE_DIRECTORY, resourceDirectory, "path")
        values.put(ConfigurationConstants.CFG_OUTPUT_BASE_DIRECTORY, context.getOutputDirectory().getAbsolutePath(), "path")
        values.put(COConstants.CVALUE_POSSIBLE_TUMOR_SAMPLE_NAME_PREFIXES, "tumor01, tumor", "string")
        values.put(COConstants.CVALUE_SAMPLE_DIRECTORY, '${inputBaseDirectory}/${pid}/${sample}/paired', "path")
        values.put(COConstants.CVALUE_SEQUENCE_DIRECTORY, '${sampleDirectory}/${run}/sequence', "path")

        [sample0, sample1].each { Sample sample ->
            def filesForSample = new COProjectsRuntimeService().loadLaneFilesForSample(context, sample)
            fileGroupSetMap[sample] = new RNAseqLaneFileGroupSetForSingleEnd(filesForSample)
        }
    }

    @Test
    void getLeftLaneFilesAsCSVs() throws Exception {
        assert fileGroupSetMap[sample0].getLeftLaneFilesAsCSVs() == "$f0l,$f1l"
    }

    @Test
    void getRightLaneFilesAsCSVs() throws Exception {
        assert fileGroupSetMap[sample0].getRightLaneFilesAsCSVs() == ""
    }

    @Test
    void getLaneFilesAlternatingWithSpaceSep() throws Exception {
        assert fileGroupSetMap[sample0].getLaneFilesAlternatingWithSpaceSep() == "$f0l $f1l"
    }

    /**
     * Will not be used in the workflow
     * @throws Exception
     */
    @Test
    @Ignore
    void getFlowCellIDsWithSpaceSepTest() throws Exception {
        assert fileGroupSetMap[sample0].getFlowCellIDsWithSpaceSep() == "BC5YE7ACXX BC6B2MACXX"
        assert fileGroupSetMap[sample1].getFlowCellIDsWithSpaceSep() == "BC5YE7ACXX BC6B2MACXX"
    }

    /**
     * This test will be ignored for single end.
     * It is not used in the RNAseqWorkflow and will fail!
     * The code in the QCPipelineScriptFileServiceHelper is not correct for single ends.
     * The problem is, that it is not easily possible to detect the lane id. For paired ends
     * the lane id is found by comparing both filenames and take everything until the first
     * differing "block" (separated by "_").
     * @throws Exception
     */
    @Test
    @Ignore
    void collectLaneIDs() throws Exception {
        assert fileGroupSetMap[sample0].getLaneIDsWithSpaceSep() == "D2826_GATCAGA_L002 D2826_GATCAGA_L002"
        assert fileGroupSetMap[sample1].getLaneIDsWithSpaceSep() == "D2826_GATCAGA_L002 D2826_GATCAGA_L002"
    }

    /**
     * Will not be used in the workflow
     * @throws Exception
     */
    @Test
    @Ignore
    void collectFlowCellAndLaneIDsWithSpaceSep() throws Exception {
        assert fileGroupSetMap[sample0].getFlowCellAndLaneIDsWithSpaceSep() == "BC5YE7ACXX D2826_GATCAGA_L002 BC6B2MACXX D2826_GATCAGA_L002"
        assert fileGroupSetMap[sample1].getFlowCellAndLaneIDsWithSpaceSep() == "BC5YE7ACXX D2826_GATCAGA_L002 BC6B2MACXX D2826_GATCAGA_L002"
    }

    @Test
    void getBamReadGroupLinesTest() throws Exception {
        def content = [
                "ID:run160319_D00133_0107_BC5YE7ACXX_D2826_GATCAGA_L002_R1_001 LB:tumor0_TEST_PID PL:ILLUMINA SM:sample_tumor0_TEST_PID PU:BC5YE7ACXX",
                "ID:run160326_D00695_0025_BC6B2MACXX_D2826_GATCAGA_L002_R1_001 LB:tumor0_TEST_PID PL:ILLUMINA SM:sample_tumor0_TEST_PID PU:BC6B2MACXX",
                "ID:run170319_D00133_0107_BC5YE7ACXX_D2826_GATCAGA_L002_R1_001 LB:tumor02_TEST_PID PL:ILLUMINA SM:sample_tumor02_TEST_PID PU:BC5YE7ACXX",
                "ID:run170326_D00695_0025_BC6B2MACXX_D2826_GATCAGA_L002_R1_001 LB:tumor02_TEST_PID PL:ILLUMINA SM:sample_tumor02_TEST_PID PU:BC6B2MACXX"
        ]

        assert fileGroupSetMap[sample0].getBamReadGroupLines() == "${content[0]} , ${content[1]}"
        assert fileGroupSetMap[sample1].getBamReadGroupLines() == "${content[2]} , ${content[3]}"
    }
}