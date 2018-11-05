package de.dkfz.b080.co.rnaseqworkflow

import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.files.COConstants
import de.dkfz.b080.co.files.Sample
import de.dkfz.roddy.RunMode
import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.config.ConfigurationConstants
import de.dkfz.roddy.core.ContextResource
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.execution.io.ExecutionService
import de.dkfz.roddy.execution.io.LocalExecutionService
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider
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
class RNAseqLaneFileGroupSetTestForPairedEnd extends ContextResource {

    @ClassRule
    public static final ContextResource contextResource = new ContextResource()

    static ExecutionContext context

    static Sample sample0 = new Sample(context, "tumor0")
    static Sample sample1 = new Sample(context, "tumor02")
    static Map<Sample, RNAseqLaneFileGroupSetForPairedEnd> fileGroupSetMap = [:]
    static boolean servicesSetup = false

    // The test directory
    static String userDir = new File(System.getProperty("user.dir")).parent
    static String pairedFolderSample0 = userDir + "/out/test/RNAseqWorkflow/resources/testdata_pairedEnd/TEST_PID/tumor0/paired/"
    String f0l = pairedFolderSample0 + "run160319_D00133_0107_BC5YE7ACXX/sequence/D2826_GATCAGA_L002_R1_001.fastq.gz"
    String f0r = pairedFolderSample0 + "run160319_D00133_0107_BC5YE7ACXX/sequence/D2826_GATCAGA_L002_R2_001.fastq.gz"
    String f1l = pairedFolderSample0 + "run160326_D00695_0025_BC6B2MACXX/sequence/D2826_GATCAGA_L002_R1_001.fastq.gz"
    String f1r = pairedFolderSample0 + "run160326_D00695_0025_BC6B2MACXX/sequence/D2826_GATCAGA_L002_R2_001.fastq.gz"

    static synchronized void setupServices() {
        if (servicesSetup) return
        ExecutionService.initializeService(LocalExecutionService.class, RunMode.CLI)
        FileSystemAccessProvider.initializeProvider(true)
        servicesSetup = true
    }

    @BeforeClass
    static void setup() {
        context = contextResource.createSimpleContext(RNAseqLaneFileGroupSetTestForPairedEnd, new Configuration(null), new COProjectsRuntimeService())

        // The setup is a bit complicated, because I want to load the fastq files from disk to simulate the full
        // fastq loading behaviour
        def resourceDirectory = LibrariesFactory.getGroovyClassLoader().getResource("resources/testdata_pairedEnd").file
        setPrivateField("inputDirectory", context, new File(resourceDirectory))

        RNAseqLaneFileGroupSetTestForPairedEnd.setupServices()

        def values = context.getConfiguration().getConfigurationValues()
        values.put(ConfigurationConstants.CFG_INPUT_BASE_DIRECTORY, resourceDirectory, "path")
        values.put(ConfigurationConstants.CFG_OUTPUT_BASE_DIRECTORY, context.getOutputDirectory().getAbsolutePath(), "path")
        values.put(COConstants.CVALUE_POSSIBLE_TUMOR_SAMPLE_NAME_PREFIXES, "tumor01, tumor", "string")
        values.put(COConstants.CVALUE_SAMPLE_DIRECTORY, '${inputBaseDirectory}/${pid}/${sample}/paired', "path")
        values.put(COConstants.CVALUE_SEQUENCE_DIRECTORY, '${sampleDirectory}/${run}/sequence', "path")

        [sample0, sample1].each { Sample sample ->
            def filesForSample = new COProjectsRuntimeService().loadLaneFilesForSample(context, sample)
            fileGroupSetMap[sample] = new RNAseqLaneFileGroupSetForPairedEnd(filesForSample)
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
     * Will not be used in the workflow
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
                "ID:run160319_D00133_0107_BC5YE7ACXX_D2826_GATCAGA_L002 LB:tumor0_TEST_PID PL:ILLUMINA SM:sample_tumor0_TEST_PID PU:BC5YE7ACXX",
                "ID:run160326_D00695_0025_BC6B2MACXX_D2826_GATCAGA_L002 LB:tumor0_TEST_PID PL:ILLUMINA SM:sample_tumor0_TEST_PID PU:BC6B2MACXX",
                "ID:run170319_D00133_0107_BC5YE7ACXX_D2826_GATCAGA_L002 LB:tumor02_TEST_PID PL:ILLUMINA SM:sample_tumor02_TEST_PID PU:BC5YE7ACXX",
                "ID:run170326_D00695_0025_BC6B2MACXX_D2826_GATCAGA_L002 LB:tumor02_TEST_PID PL:ILLUMINA SM:sample_tumor02_TEST_PID PU:BC6B2MACXX"
        ]

        assert fileGroupSetMap[sample0].getBamReadGroupLines() == "${content[0]} , ${content[1]}"
        assert fileGroupSetMap[sample1].getBamReadGroupLines() == "${content[2]} , ${content[3]}"
    }
}