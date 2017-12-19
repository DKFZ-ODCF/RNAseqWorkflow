/*
 * Copyright (c) 2017 eilslabs.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/Roddy/LICENSE.txt).
 */

package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.files.COConstants
import de.dkfz.b080.co.files.Sample
import de.dkfz.roddy.RunMode
import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.config.ConfigurationConstants
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.MockupExecutionContextBuilder
import de.dkfz.roddy.execution.io.ExecutionService
import de.dkfz.roddy.execution.io.LocalExecutionService
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider
import de.dkfz.roddy.plugins.LibrariesFactory
import groovy.transform.CompileStatic
import org.junit.BeforeClass
import org.junit.Test

import static de.dkfz.b080.co.rnaseqworkflow.Helper.setPrivateField

/**
 * Created by heinold on 05.12.16.
 */
@CompileStatic
class RNAseqLaneFileGroupSetTestForSingleCell {

    static final String SINGLECELL_SAMPLE = "singlecell"

    static ExecutionContext context = MockupExecutionContextBuilder.createSimpleContext(RNAseqLaneFileGroupSetTestForSingleCell, new Configuration(null), new COProjectsRuntimeService())
    static Sample sample = new Sample(context, SINGLECELL_SAMPLE)
    static Map<Sample, RNAseqLaneFileGroupSetForPairedEnd> fileGroupSetMap = [:]

    // The test directory
    static String userDir = new File(System.getProperty("user.dir")).parent

    static String pairedFolderSample1 = userDir + "/out/test/RNAseqWorkflow/resources/testdata/TEST_PID/${SINGLECELL_SAMPLE}/paired/"
    String g0r0 = pairedFolderSample1 + "run160319_D00133_0107_BC5YE7ACXX/sequence/D2826_GATCAGA_L002_R2_001.fastq.gz"
    String g0r1 = pairedFolderSample1 + "run160319_D00133_0107_BC5YE7ACXX/sequence/D2826_GATCAGA_L003_R2_001.fastq.gz"
    String g0r2 = pairedFolderSample1 + "run160319_D00133_0107_BC5YE7ACXX/sequence/D2826_GATCAGA_L007_R2_001.fastq.gz"
    String g0r3 = pairedFolderSample1 + "run160319_D00133_0107_BC5YE7ACXX/sequence/D2826_GATCAGA_L006_R2_001.fastq.gz"

    @BeforeClass
    static void setup() {
        // The setup is a bit complicated, because I want to load the fastq files from disk to simulate the full
        // fastq loading behaviour
        def resourceDirectory = LibrariesFactory.getGroovyClassLoader().getResource("resources/testdata").file
        setPrivateField("inputDirectory", context, new File(resourceDirectory))
        ExecutionService.initializeService(LocalExecutionService.class, RunMode.CLI)
        FileSystemAccessProvider.initializeProvider(true)
        def values = context.getConfiguration().getConfigurationValues()
        values.put(ConfigurationConstants.CFG_INPUT_BASE_DIRECTORY, resourceDirectory, "path")
        values.put(ConfigurationConstants.CFG_OUTPUT_BASE_DIRECTORY, context.getOutputDirectory().getAbsolutePath(), "path")
        values.put(COConstants.CVALUE_POSSIBLE_TUMOR_SAMPLE_NAME_PREFIXES, SINGLECELL_SAMPLE, "string")
        values.put(COConstants.CVALUE_SAMPLE_DIRECTORY, '${inputBaseDirectory}/${pid}/${sample}/paired', "path")
        values.put(COConstants.CVALUE_SEQUENCE_DIRECTORY, '${sampleDirectory}/${run}/sequence', "path")
        values.put(COConstants.FLAG_USE_SINGLE_END_PROCESSING, "true", "boolean")

        def filesForSample = new COProjectsRuntimeService().loadLaneFilesForSample(context, sample)
        fileGroupSetMap[sample] = new RNAseqLaneFileGroupSetForPairedEnd(filesForSample)
    }


    @Test
    void testSetup() {
        assert new COProjectsRuntimeService().loadLaneFilesForSample(context, sample).size() == 4
    }

    @Test
    void getSingleCellLaneFilesByRunTest() throws Exception {
        def mapByRun = fileGroupSetMap[sample].getSingleCellLaneFilesByRun()
        assert mapByRun.size() == 1
        assert mapByRun["run160319_D00133_0107_BC5YE7ACXX"].size() == 4
        assert mapByRun["run160319_D00133_0107_BC5YE7ACXX"][0].path.name == "D2826_GATCAGA_L002_R1_001.fastq.gz"
        assert mapByRun["run160319_D00133_0107_BC5YE7ACXX"][1].path.name == "D2826_GATCAGA_L002_R1_003.fastq.gz"
        assert mapByRun["run160319_D00133_0107_BC5YE7ACXX"][2].path.name == "D2826_GATCAGA_L002_R1_006.fastq.gz"
        assert mapByRun["run160319_D00133_0107_BC5YE7ACXX"][3].path.name == "D2826_GATCAGA_L002_R1_007.fastq.gz"
    }
}