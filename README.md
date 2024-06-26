## RNAseq processing workflow

Author (current): ?
Author (until 2018): Naveed Ishaque

This workflow does the primary data processing for RNAseq data: alignment, QC, read counting, reference free quantification, fusion detection. It was originally developed in the "eilslabs" at the German Cancer Research Center (DKFZ), Heidelberg. The workflow uses the workflow management system [Roddy](https://github.com/TheRoddyWMS/Roddy). Roddy manages workflow runs on cluster environments (as of 2019-12-02 this is LSF and PBS). 

### Description


#### Output

| File                                         | Description                                                                                                                                                                                                                                            |
|----------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `*_merged.mdup.bam` & `.bai`                 | Standard alignment.                                                                                                                                                                                                                                    |
| `*_chimeric_merged.mdup.bam` & `.bai`        | Chimeric alignments.                                                                                                                                                                                                                                   |
| `_chimeric_merged.junction`                  | Correspond to the file mentioned in section “6.4 Chimeric alignments in Chimeric.out.junction” in the [STAR manual](https://github.com/alexdobin/STAR/blob/2.7.10a/doc/STARmanual.pdf). The columns and their content is described in the STAR manual. |
| `_merged.mdup.bam.flagstat`                  | flag statistics of the “standard” bam file                                                                                                                                                                                                             |
| `_star_logs_and_files`                       | STAR logs                                                                                                                                                                                                                                              |
| `featureCounts/*.fpkm_tpm.featureCounts.tsv` | Gene-based expression values derived from `featureCounts`. See below.                                                                                                                                                                                  |
| `featureCounts_dexseq/*.fpkm_tpm.featureCounts.dexseq.tsv`                     | Exon-based expression values derived from `featureCounts`. See below.                                                                                                                                                                                  |
| `fusions_arriba/` | Output of the [Arriba](https://github.com/suhrig/arriba/) tool                                                                                                                                                                                         |


#### Extended Feature Counts Table

The following concerns the tables in the `featureCounts/*.fpkm_tpm.featureCounts.tsv` (gene-based) and `featureCounts_dexseq/*.fpkm_tpm.featureCounts.dexseq.tsv` (exon-based) files. The [FPKM and TPM values](https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/) are calculated in [featureCounts_2_FpkmTpm](https://github.com/DKFZ-ODCF/RNAseqWorkflow/blob/master/resources/analysisTools/rnaseqworkflow/featureCounts_2_FpkmTpm) (gene-based) and [featureCountsDexseq_2_FpkmTpm](https://github.com/DKFZ-ODCF/RNAseqWorkflow/blob/master/resources/analysisTools/rnaseqworkflow/featureCountsDexseq_2_FpkmTpm) (exon-based) from the `.s0` (unstranded), `.s1` (forward/sense), and `.s2` (reverse/antisense) output files of featureCounts.

* RPKM: Reads per kilobase gene model, per million read<br>
  <img src="https://render.githubusercontent.com/render/math?math=RPKM(g)%20=%20\frac{read\_count(g)%20/%20length(g)}{\sum_{i%20\in%20Genes}{read\_count(i)}}%20*%2010^9" style="background-color:white;">
* FPKM: Like RPKM, but for fragments with paired-end data (a pair of aligning reads is counted as one fragment)<br>
  <img src="https://render.githubusercontent.com/render/math?math=FPKM(g)%20=%20\frac{fragment\_count(g)%20/%20length(g)}{\sum_{i%20\in%20Genes}{fragment\_count(i)}}%20*%2010^9" style="background-color:white;">
* TPM: Like RPKM, but the length-weighted sum of read counts is used in the denominator<br>
  <img src="https://render.githubusercontent.com/render/math?math=TPM(gene)%20=%20\frac{RPKM(g)}{\sum_{i%20\in%20Genes}{RPKM(i)}}%20*%2010^6%20=%20\frac{read\_count(g)%20/%20length(g)}{\sum_{i%20in%20Genes}{read\_count(i)%20/%20length(i)}}%20*%2010^6" style="background-color:white;"> 


| Header     | Description |
|------------|-------------------|
| #chrom     | Chromosome contigs |
| chromStart | Start position of the gene |
| chromEnd   | End position of the gene |
| gene_id    | Ensembl gene ID. For the exon-based featureCounts file this column states all Ensembl gene IDs for which the respective exon is part of (separated by "+"). |
| score      | - |
| strand     | Strand of the gene |
| name                       | Gene name |
| exonic_length              | Length of the coding region of the gene |
| num_reads_unstranded                  | The number of reads falling into the gene region, in a non-strand specific counting mode (unstranded - `s 0` option in featureCounts).|
| num_reads_stranded               | The number(s) of reads falling into the gene region, in a strand specific counting modes. Applying the featureCounts option for standedness of `s 1`.  When using featureCounts, you can choose the `1/_stranded` or `2/_reverse_stranded` option for standedness based on the library preparation step of your stranded sequenced data. For example, the Illumina TruSeq stranded RNA kit usually produces data that fits the option `2/_reverse_stranded`. |
| num_reads_reverse_stranded               |   The number of reads falling into the gene region. The counts are based on the `s 2` option in featureCounts. To obtain the counts for the Illumina TruSeq stranded RNA kit library protocol, refer to this column as mentioned above. And also use `_reverse_stranded` columns for FPKM and TPM values from this protocol. |
| FPKM_customLibSize_unstranded                       | Unstranded FPKM calculation using the counts from the `s 0` option in feature counts. Counts for genes from rRNA, tRNA, chrX, chrY, and Mt were excluded during library size estimation.|
| FPKM_customLibSize_stranded                    | Stranded FPKM calculation based on the counts from `-s 1` option in featureCount with custom library size estimation as described above.|
| FPKM_customLibSize_reverse_stranded                    | Reverse stranded FPKM calculation based on the counts from `-s 2` option in featureCount with custom library size estimation as described above.|
| TPM_customLibSize_unstranded                        | Unstranded TPM calculation using the counts from the `s 0` option in feature counts. Counts for genes from rRNA, tRNA, chrX, chrY, and Mt were excluded during library size estimation.|
| TPM_customLibSize_stranded                     | Stranded TPM calculation was performed using the counts from the `s 1` option in feature counts with custom library size estimation as above.|
| TPM_customLibSize_reverse_stranded                     | Reverse stranded TPM calculation was performed using the counts from the `s 2` option in feature counts with custom library size estimation as above.|
| FPKM_unstranded                | Unstranded FPKM calculation based on the counts from `-s 0` option in featureCount|
| FPKM_stranded            | Stranded FPKM calculation based on the counts from `-s 1` option in featureCount|
| FPKM_reverse_stranded            | Reverse stranded FKPM calculation based on the counts from `-s 2` option in featureCount|
| TPM_unstranded                 |Unstranded TPM calculation based on the counts from `-s 0` option in featureCount |
| TPM_stranded             | Stranded TPM calculation based on the counts from `-s 1` option in featureCount|
| TPM_reverse_stranded             | Reverse stranded TPM calculation based on the counts from `-s 2` option in featureCount|


#### Example Protocol

The following is kind of a template protocol for a methods section. You will probably need to adapt it to your specific settings.

The RNAseq data were analysed with the DKFZ/ODCF RNAseq workflow (https://github.com/DKFZ-ODCF/RNAseqWorkflow, **version**; https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows, **version**; https://github.com/TheRoddyWMS/Roddy-Default-Plugin, **version**; https://github.com/TheRoddyWMS/Roddy-Base-Plugin, **version**; https://github.com/TheRoddyWMS/Roddy, **version**). The workflow performs the following analysis steps.

   NOTE: You should also list the other plugins in the dependency chain, and the Roddy version you used to ensure reproducibility. Use the short commit hash if the version you used was not tagged.

FASTQ reads for individual samples were aligned by a 2 pass alignment using the STAR aligner (**version**, https://www.ncbi.nlm.nih.gov/pubmed/23104886). Reads were aligned to a STAR index generated from the 1000 genomes assembly, gencode 19 gene models and for asjbdOverhang of 200. The alignment call parameters were -

```
--sjdbOverhang 200 \
--runThreadN 8 \
--outSAMtype BAM Unsorted SortedByCoordinate \
--limitBAMsortRAM 100000000000 \
--outBAMsortingThreadN=1 \
--outSAMstrandField intronMotif \
--outSAMunmapped Within KeepPairs \
--outFilterMultimapNmax 1 \
--outFilterMismatchNmax 5 \
--outFilterMismatchNoverLmax 0.3 \
--twopassMode Basic \
--twopass1readsN -1 \
--genomeLoad NoSharedMemory \
--chimSegmentMin 15 \
--chimScoreMin 1 \
--chimScoreJunctionNonGTAG 0 \
--chimJunctionOverhangMin 15 \
--chimSegmentReadGapMax 3 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--alignIntronMax 1100000 \
--alignMatesGapMax 1100000 \
--alignSJDBoverhangMin 3 \
--alignIntronMin 20 \
--clip3pAdapterSeq AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--readFilesCommand gunzip -c
```

Other parameters were as default, or only pertinent for particular samples (e.g. list of FASTQ files or definitions of the RG line).

Duplicate marking of the resultant main alignment file was done with sambamba (**version**, https://www.ncbi.nlm.nih.gov/pubmed/25697820) using 8 threads.

The Chimeric file was sorted using samtools sort (**version**, https://www.ncbi.nlm.nih.gov/pubmed/19505943), and then duplicates were marked using sambamba.

BAM indexes were generated using sambamba.

Quality control analysis was performed using the `samtools flagstat` command, and the rnaseqc tool (**version**, https://www.ncbi.nlm.nih.gov/pubmed/22539670) with the 1000 genomes assembly and gencode 19 gene models. Depth of Coverage analysis for rnaseqc was turned off.

Featurecounts (**version**, https://www.ncbi.nlm.nih.gov/pubmed/24227677) was used to perform gene specific read counting over exon features based on the gencode 19 gene models. Both reads of a paired fragment were used for counting and the quality threshold was set to 255 (which indicates that STAR found a unique alignment). Strand unspecific counting was used.

A custom script was used to calculate RPKM and TPM expression values. 

For total library abundance calculations, all genes on chromosomes X, Y, MT and rRNA and tRNA genes were omitted as they are likely to introduce library size estimation biases. The estimation of library complexity is done with by [rnaseqc](https://github.com/getzlab/rnaseqc/blob/master/Metrics.md), which implements the algorithm of [Picard EstimateLibraryComplexity](https://broadinstitute.github.io/picard/command-line-overview.html#EstimateLibraryComplexity).

Gene fusions were identified using the arriba algorithm (version, https://github.com/suhrig/arriba/).

### Software Versions

Usually, the software versions used for the analysis are set in one of the scripts in the  `resources/analysisTools\rnaseqworkflow\environments` directory. Which of these scripts was used for your analyses is determined by the `workflowEnvironmentScript` variable.  Dependent on that script, specific version number may also have been set in e.g. the configuration files or command-line. In this case, the actually used value can be found in the `roddyExecutionStore/exec_*/*.parameter` files. 

Currently, there are two example environment setup scripts. One of them loads a "composite" module, which makes it very hard or impossible to identify the software version just from the logging output. If such an environment is used, the only option you have is to load the "composite" module, and check the versions manually. 

The second example uses `_VERSION` variables that specify the version of the module to load. In this case, you can just do a grep on `_VERSION` on the `.parameter` file to retrieve the software versions.

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

The following is a list of citations for software used by the RNAseq workflow. Note that dependent on your configuration some tools may be unused. The actual versions used may deviate from those given in this list.

* Python 2.7.9
* [star 2.7.10a](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
* [samtools 1.6](https://www.ncbi.nlm.nih.gov/pubmed/19505943)
* [arriba 2.3.0](https://github.com/suhrig/arriba/)
* [subread 1.6.5](http://subread.sourceforge.net/) providing [featurecounts](https://www.ncbi.nlm.nih.gov/pubmed/24227677)
* [rnaseqc 1.1.8](https://www.ncbi.nlm.nih.gov/pubmed/22539670)
* [sambamba 0.6.5](https://www.ncbi.nlm.nih.gov/pubmed/25697820)
* [qualimap 2.2.1](http://qualimap.bioinfo.cipf.es/)
* [kallisto 0.46.0](https://pachterlab.github.io/kallisto/about)
* [Jemultiplexer 1.0.6](https://gbcs.embl.de/portal/tiki-index.php?page=Jemultiplexer) 

##### Conda Environment

The [Conda](https://conda.io/docs/)-environment is work in progress. The current version of the file that is delivered as an "outlook" differs from our current production environment is some aspects:

  * Jemultiplexer is yet missing from Conda
  * qualimap is not available in version 2.2.1 in Conda. Instead the environment contains version 2.2.2a.
  * Rather than R 3.0.0 the conda environment uses R 3.1.2.
  
The environment is not tested. Furthermore we noticed that Conda environments can sometimes not be restored, even if the "bioconda-legacy" channel is used. 

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
| ARRIBA_KNOWN_FUSIONS | GZipped TSV | ${hg19BaseDirectory}/tools_data/arriba/known_fusions_hg19_hs37d5_GRCh37_v2.2.1.tsv.gz | 
| ARRIBA_BLACKLIST | GZipped TSV | ${hg19BaseDirectory}/tools_data/arriba/blacklist_hg19_hs37d5_GRCh37_v2.2.1.tsv.gz |
| ARRIBA_PROTEIN_DOMAINS | GFF3 | ${hg19BaseDirectory}/tools_data/arriba/protein_domains_hg19_hs37d5_GRCh37_v2.2.1.gff3 |
| ARRIBA_CYTOBANDS | TSV | ${hg19BaseDirectory}/tools_data/arriba/cytobands_hg19_hs37d5_GRCh37_v2.2.1.tsv |

The Arriba files are available via the Arriba release tarball available at [GitHub](https://github.com/suhrig/arriba/).

#### Setting up the reference data

```bash
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz \
    && gunzip hs37d5.fa.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz -O /dev/stdout \
    | awk -F'\t' '$3!="gene"' \
    | perl -pe 's/^chr(\d+)/\1/' \
    | perl -pe 's/^chrM/MT/' \
    | perl -pe 's/^chr([XY])/\1/' \
    | perl -pe 's/ protein_id "[^"]+";//' \
    > gencode.v19.annotation_plain.noGenes.gtf

gffread \
    -w kallisto-0.43.0_1KGRef_Gencode19.noGenes.fa \
    -g hs37d5.fa \
    gencode.v19.annotation_plain.noGenes.gtf
   
kallisto-0.43.0 \
    index \
    --index=kallisto-0.43.0_1KGRef_Gencode19_k31.noGenes.index \
    kallisto-0.43.0_1KGRef_Gencode19.noGenes.fa
```

# Running the workflow

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
| RUN_KALLISTO | false | runs reference free quantification using kallisto. |
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
| GENE_MODELS_NOGENE | | For gencode < v21 this can be the same as `GENE_MODELS`. For later versions, create a file from `GENE_MODELS_NOGENE` by dropping all "gene" features. E.g. `awk '$3 != "gene"' gencode.v32.annotation_plain.gtf > gencode.v32.annotation_noGene.gtf` |
| ARRIBA_KNOWN_FUSIONS | | Only needed if arriba is used for fusion detection. ${hg19BaseDirectory}/tools_data/arriba/known_fusions_CancerGeneCensus_gencode19_2017-01-16.tsv.gz |
| ARRIBA_BLACKLIST | | Only needed if arriba is used for fusion detection. ${hg19BaseDirectory}/tools_data/arriba/blacklist_hs37d5_gencode19_2017-01-09.tsv.gz |
|--------------------|--------|-------------|

> NOTE: Please turn only at most one of the three options `RUN_KALLISTO`, `RUN_KALLISTO_FR` or `RUN_KALLISTO_RF`. These options are mutually exclusive.

## Example Call

* Change the basic stuff like `inputbase` and `outputbase` directories
set `RUN_FEATURE_COUNTS_DEXSEQ`, `RUN_RNASEQC`, `RUN_KALLISTO`, `RUN_ARRIBA`, `runFingerprinting` to `FALSE`.
* Check if your sample type is in `possibleTumorSampleNamePrefixes`
* Set your star index (`GENOME_STAR_INDEX`) and gene models (`GENE_MODELS`) parameters to what you created for (1) and (2).
* Set `GENOME_GATK_INDEX` to point to the new FASTA file (1), for which you have created a dict file. The dict file should be in the same directory as the FASTA.

## Change Log

* 4.0.0
  - major: `resources/configurationFiles/analysisRNAseq.xml` is now just a default configuration with many configuration options left blank. For your `<analysis>` tags in your project XMLs the analysis names can be changed to use the following defaults:
     * GRCh37, mm10: [RNAseqAnalysis](resources/configurationFiles/analysisRNAseq.xml)
     * GRCh38-specific: [RNAseqAnalysisGRCh38](resources/configurationFiles/analysisRNAseqGRCh38.xml)
  - major: Update default software versions for all assemblies
    * Arriba 2.3.0
    * STAR 2.7.10a (from 2.5.3a for hg37/mm10 and 2.7.6a for hg38)
    * Kallisto 0.46.0 (from 0.42.0)
    * Samtools 1.9 (from 1.6)
    * HTSlib 1.9 (from 1.6)
    * subread 1.6.5 (from 1.5.1): The previous version produces occasional segmentation faults (related to extreme optimization option `-O6`) but otherwise produces the same results. Both versions produced exactly the same `featureCounts` in multiple tests.
  - major: Gene model gencode version 39/43 was used to create the STAR and Kallisto indexes
  - major: Update default paths to new ngs_share at ODCF.
  - major: GRCh38 STAR and Kallisto indexes are based on `refmake` workflow.
  - major: Updated the Conda `environment.yaml`.
    * Note for users at the DKFZ cluster: The software versions used in Conda environment do not exactly match the versions in the default configuration used for the cluster's module system.
  - major: Column rename in featureCounts table:
    - "FPKM_no_mt_rrna_trna_chrxy{,_fw,_rv}" -> "FPKM_customLibSize{_unstranded, _stranded, _reverse_stranded}"
    - "TPM_no_mt_rrna_trna_chrxy{,_fw,_rv}" -> "TPM_customLibSize{_unstranded, _stranded, _reverse_stranded}"
    - "FPKM_standard{,_fw,_rv}" -> "FPKM{_unstranded, _stranded, _reverse_stranded}"
    - "TPM_standard{,_fw,_rv}" -> "TPM{_unstranded, _stranded, _reverse_stranded}"

* 3.0.0 [19th Oct 2021]
  - major: Column rename in feature counts table:
    - "FPKM{,_fw,_rv}" -> "FPKM_no_mt_rrna_trna_chrxy{,_fw,_rv}"
    - "TPM{,_fw,_rv}" -> "TPM_no_mt_rrna_trna_chrxy{,_fw,_rv}"
    - "FPKM_legacy{,_fw,_rv}" -> "FPKM_standard{,_fw,_rv}"
    - "TPM_legacy{,_fw,_rv}" -> "TPM_standard{,_fw,_rv}"

* 2.1.0 (2.0.3-deprecated) [26th Mar 2021]
  - minor: Added GRCh38 genome support; version changes only affect hg38
    - Reference genome (core_ref_GRCh38_hla_decoy_ebv.tar.gz) was downloaded from ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/ and Illumina PhiX genome was added.
    - Added GRCh38 specific configs to a `resources/configurationFiles/analysisRNAseqGRCh38.xml`
    - minor: STAR to 2.7.6a
    - minor: Updated Kallisto to 0.46.0
    - minor: Gene model gencode version 31 was used to create the STAR and Kallisto indexes
  - minor: Update Samtools to 1.9
  - minor: Update HTSlib to 1.9
  - 2.0.3 was deprecated to correct the incorrect version number. 

* 2.0.2 [27th Jan 2021]
  - patch: Rnaseqc requires "transcript_id" attributes also for "gene" features. Gencode versions >= v21, however, don't have a "transcript_id" attribute in the "gene" features. You can now provide a Gencode annotation without gene features via the `GENE_MODELS_NOGENE` variable.
  - patch: Renamed `HIPO2_rnaseq_processing.sh` to `rnaseq_processing.sh`.

* 2.0.1 [16th Jul 2020]
  - patch: bugfix: Variable `SEQUENCING_PROTOCOL` was statically overwritten to be "paired". In a mode where the plugin retrieves the files from the filesystem (i.e. without `FASTQ_LIST`), the workflow therefore searched for FASTQs below the "paired" rather than "single" directory.

* 2.0.0 [13th Mar 2020]
  - major: Update Arriba to version 1.2.0
  - patch: Update Subread (featureCounts) to 1.5.3
  - patch: Added draft Conda environment (incomplete)

* 1.3.0-2 [26th Nov 2019]
  - Removed the single-quotes around `${ADAPTER_SEQ}` in `--clip3pAdapterSeq` again. STAR uses non-standard way of parsing parameters and manages to get all adapters. With quotes the adapters get also, and it is unclear what STAR does with them, except that it does not complain about a configuration error and that it also does not complain with even more severe misconfigurations, such as other non-DNA sequences as adapter sequences. The manual also does not use quoted parameter arguments, so no-quotes is conform to this STAR-specific CLI parameter handling pattern.

* 1.3.0-1 [7th Nov 7 2018]
  - patch: Added single quotes around `$ADAPTER_SEQ` parameter in `--clip3pAdapterSeq` to allow for separate first and second read adapters.
  - patch: bugfix: Changed all plain "$varName" references in the configuration file to "${varName}" references. Roddy only knows the latter.

* 1.3.0 [9th Apr 2018]
  - minor: works with Roddy 3

* 1.2.23-2 [15th Feb 2018]
  - Modified software defaults to samtools 1.6, star 2.5.2b -> 2.5.3a, arriba 0.8 -> 0.12.
  - Added exception for loading htslib if samtools version is too low. 
  - Modified output file check of BAM file
  
* 1.2.23-1 [26th Jan 2018]
  - Changed default samtools from v1.3.1 to v1.6

* 1.2.23 [25th Jan 2018]
  - Update to Roddy 3.5 that fixes a non-quoting error with the Bash environment and spaces in referenced variables (relevant for adapter variables).
    Note that if you did set your adapter sequences by reference to another variable, e.g. "ADAPTER_SEQ=${ADAPTER_SEQ_TRUSEQ_LT_HT}", the previous
    Roddy versions inserted the second adapter correctly, but did not quote it. Therefore, Bash interpreted the `ADAPTER_SEQ` variable as array, and 
    failed -- for at least Bash < 4.2 -- to export both adapters correctly to the called job script.

* 1.2.22-8 [26th Nov 2019]
  - Removed the single-quotes around `${ADAPTER_SEQ}` in `--clip3pAdapterSeq` again. STAR uses non-standard way of parsing parameters and manages to get all adapters. With quotes the adapters get also quoted, and it is unclear what STAR does with them, except that it does not complain about a configuration error and that it also does not complain with even more severe misconfigurations, such as other non-DNA sequences as adapter sequences. The manual also does not use quoted parameter arguments, so no-quotes is conform to this STAR-specific CLI parameter handling pattern.

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
  - Added variable 
  
  
  _COUNT_CORES to bypass segfault with too many core
  
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

