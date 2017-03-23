package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.files.COConstants
import de.dkfz.b080.co.files.Sample
import de.dkfz.b080.co.qcworkflow.QCPipeline
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

import java.lang.reflect.Field

/**
 * Created by heinold on 05.12.16.
 */
@CompileStatic
class RNAseqLaneFileGroupSetTest {

    static ExecutionContext context = MockupExecutionContextBuilder.createSimpleContext(RNAseqLaneFileGroupSetTest, new Configuration(null), new COProjectsRuntimeService())
    static Sample sample0 = new Sample(context, "tumor0")
    static Sample sample1 = new Sample(context, "tumor02")
    static Map<Sample, RNAseqLaneFileGroupSet> fileGroupSetMap = [:]

    // The test directory
    static String userDir = new File(System.getProperty("user.dir")).parent
    static String pairedFolderSample0 = userDir + "/out/test/RNAseqWorkflow/resources/testdata/TEST_PID/tumor0/paired/"
    String f0l = pairedFolderSample0 + "run160319_D00133_0107_BC5YE7ACXX/sequence/D2826_GATCAGA_L002_R1_001.fastq.gz"
    String f0r = pairedFolderSample0 + "run160319_D00133_0107_BC5YE7ACXX/sequence/D2826_GATCAGA_L002_R2_001.fastq.gz"
    String f1l = pairedFolderSample0 + "run160326_D00695_0025_BC6B2MACXX/sequence/D2826_GATCAGA_L002_R1_001.fastq.gz"
    String f1r = pairedFolderSample0 + "run160326_D00695_0025_BC6B2MACXX/sequence/D2826_GATCAGA_L002_R2_001.fastq.gz"

    private static setPrivateField(String name, Object object, Object value) {
        Field f = null
        Class cls = object.class
        while (!f && cls) {
            try {
                f = cls.getDeclaredField(name)
            } catch (Exception ex) {
            }
            cls = cls.superclass
        }
        assert f
        f.setAccessible(true)
        f.set(object, value)
    }

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
        values.put(COConstants.CVALUE_POSSIBLE_TUMOR_SAMPLE_NAME_PREFIXES, "tumor01, tumor", "string")
        values.put(COConstants.CVALUE_SAMPLE_DIRECTORY, '${inputBaseDirectory}/${pid}/${sample}/paired', "path")
        values.put(COConstants.CVALUE_SEQUENCE_DIRECTORY, '${sampleDirectory}/${run}/sequence', "path")

        [sample0, sample1].each { Sample sample ->
            def filesForSample = new COProjectsRuntimeService().loadLaneFilesForSample(context, sample)
            fileGroupSetMap[sample] = new RNAseqLaneFileGroupSet(filesForSample)
        }
    }

    @Test
    void getLeftLaneFilesAsCSVs() throws Exception {
        assert fileGroupSetMap[sample0].getLeftLaneFilesAsCSVs() == "$f0l,$f1l"
    }

    @Test
    void getRightLaneFilesAsCSVs() throws Exception {
        assert fileGroupSetMap[sample0].getRightLaneFilesAsCSVs() == "$f0r,$f1r"
    }

    @Test
    void getLaneFilesAlternatingWithSpaceSep() throws Exception {
        assert fileGroupSetMap[sample0].getLaneFilesAlternatingWithSpaceSep() == "$f0l $f0r $f1l $f1r"
    }

    @Test
    void getFlowCellIDsWithSpaceSepTest() throws Exception {
        assert fileGroupSetMap[sample0].getFlowCellIDsWithSpaceSep() == "BC5YE7ACXX BC6B2MACXX"
        assert fileGroupSetMap[sample1].getFlowCellIDsWithSpaceSep() == "BC5YE7ACXX BC6B2MACXX"
    }

    @Test
    void collectLaneIDs() throws Exception {
        assert fileGroupSetMap[sample0].getLaneIDsWithSpaceSep() == "D2826_GATCAGA_L002 D2826_GATCAGA_L002"
        assert fileGroupSetMap[sample1].getLaneIDsWithSpaceSep() == "D2826_GATCAGA_L002 D2826_GATCAGA_L002"
    }

    @Test
    void collectFlowCellAndLaneIDsWithSpaceSep() throws Exception {
        assert fileGroupSetMap[sample0].getFlowCellAndLaneIDsWithSpaceSep() == "BC5YE7ACXX D2826_GATCAGA_L002 BC6B2MACXX D2826_GATCAGA_L002"
        assert fileGroupSetMap[sample1].getFlowCellAndLaneIDsWithSpaceSep() == "BC5YE7ACXX D2826_GATCAGA_L002 BC6B2MACXX D2826_GATCAGA_L002"
    }

    @Test
    void getBamReadGroupLinesTest() throws Exception {
        def content = [
                "ID:run160319_D00133_0107_BC5YE7ACXX_D2826_GATCAGA_L002 LB:tumor0_TEST_PID PL:ILLUMINA SM:sample_tumor0_TEST_PID PU:BC5YE7ACXX",
                "ID:run160326_D00695_0025_BC6B2MACXX_D2826_GATCAGA_L002 LB:tumor0_TEST_PID PL:ILLUMINA SM:sample_tumor0_TEST_PID PU:BC6B2MACXX",
                "ID:run170319_D00133_0107_BC5YE7ACXX_D2826_GATCAGA_L002 LB:tumor02_TEST_PID PL:ILLUMINA SM:sample_tumor02_TEST_PID PU:BC5YE7ACXX",
                "ID:run170326_D00695_0025_BC6B2MACXX_D2826_GATCAGA_L002 LB:tumor02_TEST_PID PL:ILLUMINA SM:sample_tumor02_TEST_PID PU:BC6B2MACXX"
        ]

        assert fileGroupSetMap[sample0].getBamReadGroupLines() == "${content[0]} , ${content[1]}"
        assert fileGroupSetMap[sample1].getBamReadGroupLines() == "${content[2]} , ${content[3]}"
    }
}