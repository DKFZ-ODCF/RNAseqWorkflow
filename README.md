## RNAseq processing workflow

Author (current): ?
Author (until 2018): Naveed Ishaque

This workflow does the primary data processing for RNAseq data: alignment, QC, read counting, reference free quantification, fusion detection. It was originally developed in the "eilslabs" at the German Cancer Research Center (DKFZ), Heidelberg. The workflow uses the workflow management system [Roddy](https://github.com/TheRoddyWMS/Roddy). Roddy manages workflow runs on cluster environments (as of 2019-12-02 this is LSF and PBS). 

### Description

The following is kind of a template protocol for a methods section. You will probably need to adapt it to your specific settings and we can not guarantee to always keep this up to date.

#### Example Protocol

The RNAseq data were analysed with the DKFZ/ODCF RNAseq workflow (https://github.com/DKFZ-ODCF/RNAseqWorkflow, version; https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows, version; https://github.com/TheRoddyWMS/Roddy-Default-Plugin, version; https://github.com/TheRoddyWMS/Roddy-Base-Plugin, version; https://github.com/TheRoddyWMS/Roddy, version). The workflow performs the following analysis steps.

   NOTE: You should also list the other plugins in the dependency chain and the Roddy version you used to ensure reproducibility. Use the short commit hash if the version you used was not tagged.

FASTQ reads for individual samples were aligned by a 2 pass alignment using the STAR aligner (version, https://www.ncbi.nlm.nih.gov/pubmed/23104886). Reads were aligned to a STAR index generated from the 1000 genomes assembly, gencode 19 gene models and for asjbdOverhang of 200. The alignment call parameters were -

```
--sjdbOverhang 200 --runThreadN 8 --outSAMtype BAM Unsorted SortedByCoordinate --limitBAMsortRAM 100000000000 --outBAMsortingThreadN=1 --outSAMstrandField intronMotif --outSAMunmapped Within KeepPairs --outFilterMultimapNmax 1 --outFilterMismatchNmax 5 --outFilterMismatchNoverLmax 0.3 --twopassMode Basic --twopass1readsN -1 --genomeLoad NoSharedMemory --chimSegmentMin 15 --chimScoreMin 1 --chimScoreJunctionNonGTAG 0 --chimJunctionOverhangMin 15 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 --alignIntronMax 1100000 --alignMatesGapMax 1100000 --alignSJDBoverhangMin 3 --alignIntronMin 20 --clip3pAdapterSeq AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --readFilesCommand gunzip -c
```

Other parameters were as default, or only pertinent for particular samples (e.g. list of FASTQ files or definitions of the RG line).

Duplicate marking of the resultant main alignment file was performed using sambamba (version, https://www.ncbi.nlm.nih.gov/pubmed/25697820) using 8 threads.

The Chimeric file was sorted using samtools sort (version, https://www.ncbi.nlm.nih.gov/pubmed/19505943), and then duplicates were marked using sambamba.

BAM indexes were generated using sambamba.

Quality control analysis was performed using the `samtools flagstat` command, and the rnaseqc tool (version, https://www.ncbi.nlm.nih.gov/pubmed/22539670) with the 1000 genomes assembly and gencode 19 gene models. Depth of Coverage analysis for rnaseqc was turned off.

Featurecounts (version, https://www.ncbi.nlm.nih.gov/pubmed/24227677​) was used to perform gene specific read counting over exon features based on the gencode 19 gene models. Both reads of a paired fragment were used for counting and the quality threshold was set to 255 (which indicates that STAR found a unique alignment). Strand unspecific counting was used.

A custom script was used to calculate RPKM and TPM expression values. For total library abundance calculations, all genes on chromosomes X, Y, MT and rRNA and tRNA genes were omitted as they are likely to introduce library size estimation biases.

Gene fusions were identified using the arriba algorithm (version, https://github.com/suhrig/arriba/).

### Installation of ...

#### ... the workflow management system

Please see the [Roddy documentation](https://roddy-documentation.readthedocs.io/) for details on the installation of Roddy and the [DefaultPlugin](https://github.com/TheRoddyWMS/Roddy-Default-Plugin) and [PluginBase](https://github.com/TheRoddyWMS/Roddy-Base-Plugin) plugins.

#### ... plugins

The RNA-seq plugin depends on the [AlignmentAndQCWorkflows](https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows) plugin and its dependencies. You need to install these plugins manually. Although that is somewhat inconvenient, it is also not terribly hard. The key is to understand that each plugin contains a `buildinfo.txt` file in its root directory that lists the compatible major and the minimum compatible minor versions of all directly dependent plugins. After having installed the direct dependencies, one needs to check the indirect dependencies, and so forth, until the [DefaultPlugin](https://github.com/TheRoddyWMS/Roddy-Default-Plugin) and [PluginBase](https://github.com/TheRoddyWMS/Roddy-Base-Plugin) are installed. The Roddy version to use is the same major version number and at least highest minor version number of all plugins in the chain. 

#### ... the workflow

You can install the workflow by downloading the release tarball from Github Releases, or by cloning the repository. In both cases the workflow plugin root directory should be in a direct subdirectory of your plugins directory that you register, e.g. via the `applicationProperties.ini` of Roddy.

Create a plugin directory, e.g. called `plugins`, and clone the repository into that directory. The following commands create the plugins directory and clone the repository with checked out branch "ReleaseBranch_1.3.0-2" into a correspondingly named subdirectory:

```bash
mkdir plugins
cd plugins
git clone -- https://https://github.com/DKFZ-ODCF/RNAseqWorkflow RNAseqWorkflow_2.0.0
git -C RNAseqWorkflow_2.0.0 checkout 2.0.0
```

Note that the workflow directory can be suffixed by the version tag to allow for the installation of multiple versions of the workflow. The workflow additionally needs the [AlignmentAndQCWorkflow](https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows) as the next plugin in the dependency chain.

#### ... the software stack

##### Manual Installation

Software stack version and citations for RNAseq workflow (not all tools may be used, dependent on your configuration):

* Python 2.7.9
* [star 2.5.3a](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
* [samtools 1.6](https://www.ncbi.nlm.nih.gov/pubmed/19505943)
* [arriba 1.2.0](https://github.com/suhrig/arriba/)
* [subread 1.5.3](http://subread.sourceforge.net/) providing [featurecounts 1.5.3](https://www.ncbi.nlm.nih.gov/pubmed/24227677​)
* [rnaseqc 1.1.8](https://www.ncbi.nlm.nih.gov/pubmed/22539670)
* [sambamba 0.6.5](https://www.ncbi.nlm.nih.gov/pubmed/25697820)
* [qualimap 2.2.1](http://qualimap.bioinfo.cipf.es/)
* [kallisto 0.43.0](https://pachterlab.github.io/kallisto/about)
* [Jemultiplexer 1.0.6](https://gbcs.embl.de/portal/tiki-index.php?page=Jemultiplexer) 

##### Conda Environment

The [Conda](https://conda.io/docs/)-environment is work in progress. The current version of the file than is delivered as an outlook differs from our current production environment is some aspects:

  * Jemultiplexer is yet missing from Conda
  * qualimap is not available in version 2.2.1 in Conda. Instead the environment contains version 2.2.2a.
  * Rather than R 3.0.0 the conda environment uses R 3.1.2.
  
The environment is not tested. 

#### ... the reference data

The reference data is configurable via the following configuration values. 

| Configuration value | Description | Default Path |
|---------------------|-------------|---------| 
| GENOME_FA |  FASTA assembly file | ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz |
| GENOME_GATK_INDEX | The path to the GATK index directory and name. The value may end in `.fa` just because this may be the prefix of index | |
| GENE_MODELS | GTF with gene models. | https://www.encodeproject.org/files/gencode.v19.annotation/@@download/gencode.v19.annotation.gtf.gz |
| GENE_MODELS_EXCLUDE | GTF | ${databaseDirectory}/gencode/gencode19/gencode.v19.annotation_plain.chrXYMT.rRNA.tRNA.gtf |
| GENE_MODELS_DEXSEQ | GFF | ${databaseDirectory}/gencode/gencode19/gencode.v19.annotation_plain.dexseq.gff |
| GENE_MODELS_GC | | ${databaseDirectory}/gencode/gencode19/gencode.v19.annotation_plain.transcripts.autosomal_transcriptTypeProteinCoding_nonPseudo.1KGRef.gc |
| GENOME_STAR_INDEX_50 | Path to STAR index for read length 50 pb | ${indexDirectory}/STAR/STAR_2.5.2b_1KGRef_PhiX_Gencode19_50bp |
| GENOME_STAR_INDEX_100 | Path to STAR index for read length 100 pb | ${indexDirectory}/STAR/STAR_2.5.2b_1KGRef_PhiX_Gencode19_100bp  |
| GENOME_STAR_INDEX_200 | Path to STAR index for read length 200 pb | ${indexDirectory}/STAR/STAR_2.5.2b_1KGRef_PhiX_Gencode19_200bp |
| GENOME_KALLISTO_INDEX | |  ${indexDirectory}/kallisto/kallisto-0.43.0_1KGRef_Gencode19_k31/kallisto-0.43.0_1KGRef_Gencode19_k31.noGenes.index |
| ARRIBA_KNOWN_FUSIONS | GZipped TSV | ${hg19BaseDirectory}/tools_data/arriba/known_fusions_CancerGeneCensus_gencode19_2017-05-11.tsv.gz | 
| ARRIBA_BLACKLIST | GZipped TSV | ${hg19BaseDirectory}/tools_data/arriba/blacklist_hg19_hs37d5_GRCh37_2018-11-04.tsv.gz |
| ARRIBA_PROTEIN_DOMAINS | GFF3 | ${hg19BaseDirectory}/tools_data/arriba/protein_domains_hg19_hs37d5_GRCh37_2019-07-05.gff3 |
| ARRIBA_CYTOBANDS | TSV | ${hg19BaseDirectory}/tools_data/arriba/cytobands_hg19_hs37d5_GRCh37_2018-02-23.tsv |

The Arriba files are available via the Arriba release tarball available at [GitHub](https://github.com/suhrig/arriba/).



### Running the workflow

## Configuration Values

The authoritative documentation of the configuration is the configuration itself at [resources/configurationFiles/analysisRNAseq.xml](https://www.github.com/RNAseqWorkflow/master/resources/configurationFiles/analysisRNAseq.xml). Please look at the specific version of the workflow that you are using! 

The following is merely on overview over the most important parameters. 

| Configuration value | Default | Description |
|---------------------|---------|-------------| 
| RUN_STAR	| true | runs the primary alignment using star |
| RUN_FEATURE_COUNTS | true | runs the transcript based read counting using feature counts |
| RUN_FEATURE_COUNTS_DEXSEQ	| false | runs the exon based read counting using feature counts |
| RUN_RNASEQC | true | runs RNAseQC for quality metric |
| RUN_QUALIMAP | false | runs QualiMap2 QC. This is NOT recommended as it uses all cores on the node |
| RUN_KALLISTO | false | runs reference free quantification using kallisto |
| RUN_KALLISTO_FR | false | runs reference free quantification using kallisto, forward strand specific |
| RUN_KALLISTO_RF | false | runs reference free quantification using kallisto, reverse strand specific |
| RUN_ARRIBA | true | run fusion detection using arriba |
| RUN_QCJSON | true | creates a json file with QC metrics from RNAseQC and flagstats |
| RUN_CLEANUP| true | cleans up intermediate files |
| disableDoC_GATK | false | disables Depth of Coverage calculations in RNAseQC. This is required for some x-ten samples (and I don't know why) |
| LOAD_MODULE | true | use software stack as defined by modules |
| MODULE_ENV | HIPO2_rna/v1 | Environment module to load; only applies to DKFZ/ODCF cluster environment |
| TEST_RUN | false | perform a test run without producing any ouput |
| DO_FIRST | "" | a place to enter some commands to run before the workflow begins |
| ADAPTER_SEQ | ${ADAPTER_SEQ_TRUSEQ_LT_HT} | Sequences of the adapters to trim off by STAR. This is usually one or two (space delimited) DNA-sequences. |
| GENOME_* | | Multiple configuration values defining the paths to STAR indices and genome reference files. Please see the variables starting with "GENOME_" in the configuration to get an idea of the default directory structure and files. The layout defined in the default in the plugin XML can be specified using the `indexDirectory` variable. |
| GENE_MODELS | | By default the gencode v19 annotations. See the XML for further information | 
| ARRIBA_KNOWN_FUSIONS | | Only needed if arriba is used for fusion detection. ${hg19BaseDirectory}/tools_data/arriba/known_fusions_CancerGeneCensus_gencode19_2017-01-16.tsv.gz |
| ARRIBA_BLACKLIST | | Only needed if arriba is used for fusion detection. ${hg19BaseDirectory}/tools_data/arriba/blacklist_hs37d5_gencode19_2017-01-09.tsv.gz |
|--------------------|--------|-------------|

## Example Call

* Change the basic stuff like inputbase and outputbase directories
set RUN_FEATURE_COUNTS_DEXSEQ, RUN_RNASEQC, RUN_KALLISTO, RUN_ARRIBA, runFingerprinting to FALSE.
* Check if your sample type is in possibleTumorSampleNamePrefixes
G Set your star index (GENOME_STAR_INDEX) and gene models (GENE_MODELS) parameters to what your created for (1) and (2).
* Set GENOME_GATK_INDEX to point to the new FASTA file (1), for which you have created a dict file. The dict file should be in the same directory as the FASTA.

## Change Log

* 2.0.0 [Feb 2020]
  - Update Arriba to version 1.2.0
  - Added draft Conda environment (incomplete)

* 1.3.0-2 [26th Nov 2019]
  - Removed the single-quotes around `${ADAPTER_SEQ}` in `--clip3pAdapterSeq` again. STAR uses non-standard way of parsing parameters and manages to get all adapters. With quotes the adapters get also quoted and it is unclear what STAR does with them, except that it does not complain about a configuration error and that it also does not complain with even more severe misconfigurations, such as other non-DNA sequences as adapter sequences. The manual also does not use quoted parameter arguments, so no-quotes is conform to this STAR-specific CLI parameter handling pattern.

* 1.3.0-1 [7th Nov 7 2018]
  - Added single quotes around `$ADAPTER_SEQ` parameter in `--clip3pAdapterSeq` to allow for separate first and second read adapters.
  - Changed all plain "$varName" references in the configuration file to "${varName}" references. Roddy only knows the latter.

* 1.3.0 [9th Ap 2018]
  - works with Roddy 3

* 1.2.23-2 [15th Feb 2018]
  - Modified software defaults to samtools 1.6, star 2.5.3a, arriba 0.12.
  - Added exception for loading htslib if samtools version is too low. 
  - Modified output file check of BAM file
  
* 1.2.23-1 [26th Jan 2018]
  - changed default samtools to v1.6 from v1.3.1

* 1.2.23 [25th Jan 2018]
  - Update to Roddy 3.5 that fixes a non-quoting error with the Bash environment and spaces in referenced variables (relevant for adapter variables).
    Note that if you did set your adapter sequences by reference to another variable, e.g. "ADAPTER_SEQ=${ADAPTER_SEQ_TRUSEQ_LT_HT}", the previous
    Roddy versions inserted the second adapter correctly, but did not quote it. Therefore, Bash interpreted the `ADAPTER_SEQ` variable as array, and 
    failed -- for at least Bash < 4.2 -- to export both adapters correctly to the called job script.

* 1.2.22-8 [26th Nov 2019]
  - Removed the single-quotes around `${ADAPTER_SEQ}` in `--clip3pAdapterSeq` again. STAR uses non-standard way of parsing parameters and manages to get all adapters. With quotes the adapters get also quoted and it is unclear what STAR does with them, except that it does not complain about a configuration error and that it also does not complain with even more severe misconfigurations, such as other non-DNA sequences as adapter sequences. The manual also does not use quoted parameter arguments, so no-quotes is conform to this STAR-specific CLI parameter handling pattern.

* 1.2.22-7 [7th Nov 2018]
  - Fixed a number of variable references that were without braces (now everywhere `${...}` is used to allow Roddy to resolve the references and order the parameter file correctly)
  - Updated tests to use ContextResource
  - Added single quotes around `$ADAPTER_SEQ` parameter in `--clip3pAdapterSeq` to allow for separate first and second read adapters.

* 1.2.22-6 [21st Feb 2018]
  - Lifted to Roddy 3 (stable)

* 1.2.22-5 [1st Feb 2018]
  - Added R (3.0.0) module loading to environment script

* 1.2.22-4 [31st Jan 2018]
  - Output-filename variable name change 

* 1.2.22-3 [26th Jan 2018]
  - Samtools version update from 1.3.1 to 1.6

* 1.2.22-1 [25th Jan 2018]
  - Lift to Roddy 2.4 (development version!)
  - Added module loading environment scripts for DKFZ/ODCF LSF and PBS clusters

* 1.0.22-2 [21st Sep 2017]
  - Modified default scratch to avoid crosstalk between multi sample PIDs

* 1.0.22-1 [15th Sep 2017]
  - Added variable FEATURE_COUNT_CORES to bypass segfault with too many core
  
* 1.0.22 [21st Aug 2017]
  - Fixed errors is SE qc json file (removed empty value and NaN values)

* 1.0.21 [21st Aug 2017]
  - Added single end support via "useSingleEndProcessing:true"
  - featureCounts has issues with the native bam file produced by STAR. We now use the mdup bam file and read sort the PE file (SE are read as is)
  - readgroups generated for the single end FASTQ files, which are similar to the PE readgroups
  - single end mode has 3 fewer QC values in the QC JSON. The values for "properlyPaired", "properlyPairedPercentage", "singletons", "singletonsPercentage" are no longer reported for SE mode

* 1.0.11 [12th Jun 2017]
  - Added README
  - Added tmp-directory parameter to `sambamba markdup` (it previously wrote to /temp, but now it write to `$SCRATCH`)
  - Changed default scratch to `outputAnalysisDir`
 
* 1.0.10 [9th Jun 2017]
  - Added `disableDoC_GATK` flag for some samples where DoC calculations fail in RNAseQC
  - Added additional flags for runnings kallisto in strand specific mode
  - Corrected genome index for kallisto (previously it used gene and transcript definitions from the GTF file, now it correctly only uses transcripts)

* 1.0.8 [30th Mar 2017]
  - First official release

