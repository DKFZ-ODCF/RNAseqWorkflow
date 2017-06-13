== Description

RNAseq processing workflow
Author: Naveed Ishaque

This workflows does the primary data procing for RNAseq data: alignment, QC, read counting, reference free quantification, fusion detection.

Currently Jeongbin Park is implementing a single cell analysis workflow based on this one.

SOFTWARE
```
PYTHON_BINARY   /ibios/tbi_cluster/13.1/x86_64/bin/python
STAR_BINARY     /ibios/tbi_cluster/13.1/x86_64/STAR/STAR-2.5.2b/bin/STAR
FEATURECOUNTS_BINARY    /ibios/tbi_cluster/13.1/x86_64/subread/1.5.1/bin/featureCounts
SAMBAMBA_BINARY /ibios/tbi_cluster/13.1/x86_64/sambamba/sambamba-0.6.5/sambamba
SAMTOOLS_BINARY /ibios/tbi_cluster/13.1/x86_64/bin/samtools-1.3.1
RNASEQC_BINARY  /ibios/tbi_cluster/13.1/x86_64/RNA-SeQC/RNA-SeQC-1.1.8/rnaseqc
KALLISTO_BINARY /ibios/tbi_cluster/13.1/x86_64/bin/kallisto-0.43.0
QUALIMAP_BINARY /ibios/tbi_cluster/13.1/x86_64/qualimap/qualimap-2.2.1/qualimap
ARRIBA_BINARY   /ibios/tbi_cluster/13.1/x86_64/arriba/arriba-0.8/bin/arriba
ARRIBA_READTHROUGH_BINARY       /ibios/tbi_cluster/13.1/x86_64/arriba/arriba-0.8/bin/extract_read-through_fusions
```
MODULES
``
module load HIPO2_rna/v1

module-whatis    HIPO2/v1 - This module will load further modules defined for HIPO2 in specific versions, namely: star-2.5.2b, featureCounts(subread-1.4.1), sambamba-0.6.5, rnaseqc-1.1.8, kallisto-0.43.0, qualimap-2.2.1, arriba-0.8 and 
samtools-1.3.1 
module           load STAR/2.5.2b 
module           load subread/1.5.1 
module           load sambamba/0.6.5 
module           load rnaseqc/1.1.8 
module           load kallisto/0.43.0 
module           load qualimap/2.2.1 
module           load arriba/0.8 
module           load samtools/1.3.1 
module           load Jemultiplexer/1.0.6 
module           load python/2.7.9 
``
== Run flags / switches / passable values
``
RUN_STAR			runs the primary alignment using star [true]
RUN_FEATURE_COUNTS		runs the transcript based read counting using feature counts [true]
RUN_FEATURE_COUNTS_DEXSEQ	runs the exon based read counting using feature counts [false]
RUN_RNASEQC			runs RNAseQC for quality metric [true]
RUN_QUALIMAP			runs QualiMap2 QC. This is NOT recommended as it uses all cores on the node [false]
RUN_KALLISTO			runs reference free quantification using kallisto [false]
RUN_KALLISTO_FR			runs reference free quantification using kallisto, forward strand specific [false]
RUN_KALLISTO_RF			runs reference free quantification using kallisto, reverse strand specific [false]
RUN_ARRIBA 			run fusion detection using arriba [true]
RUN_QCJSON			creates a json file with QC metrics from RNAseQC and flagstats [true]
RUN_CLEANUP			cleans up intermediate files [true]

disableDoC_GATK			disabbles Depth of Coverage calculations in RNAseQC. This is required for some x-ten samples (and I don't know why) [false]

LOAD_MODULE			use software stack as defined by modules [true]
MODULE_ENV			which module to load [HIPO2_rna/v1]
TEST_RUN			perform a test run without producing any ouput [false]

DO_FIRST			a place to enter some commands to run before the workflow begins [""]
``
== Adapter trimming params (for star)
``
ADAPTER_SEQ_TRUSEQ_LT_HT	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
ADAPTER_SEQ_TRUSEQ_DNAME	AGATCGGAAGAGCACACGTCTGAAC AGATCGGAAGAGCGTCGTGTAGGGA
ADAPTER_SEQ_TRUSEQ_SRNA		TGGAATTCTCGGGTGCCAAGG
ADAPTER_SEQ_TRUSEQ_RIBO		AGATCGGAAGAGCACACGTCT
ADAPTER_SEQ_NEXTERA		CTGTCTCTTATACACATCT
ADAPTER_SEQ_NEXTERA_MP		CTGTCTCTTATACACATCT AGATGTGTATAAGAGACAG
ADAPTER_SEQ			${ADAPTER_SEQ_TRUSEQ_LT_HT}
``
== Reference files
``
GENOME_FA			${indexDirectory}/bwa/bwa06_1KGRef_Phix/hs37d5_PhiX.fa
GENOME_GATK_INDEX		${indexDirectory}/bwa/bwa06_1KGRef_Phix/hs37d5_PhiX.fa
GENOME_STAR_INDEX_50		${indexDirectory}/STAR/STAR_2.5.2b_1KGRef_PhiX_Gencode19_50bp
GENOME_STAR_INDEX_100		${indexDirectory}/STAR/STAR_2.5.2b_1KGRef_PhiX_Gencode19_100bp
GENOME_STAR_INDEX_200		${indexDirectory}/STAR/STAR_2.5.2b_1KGRef_PhiX_Gencode19_200bp
GENOME_KALLISTO_INDEX		${indexDirectory}/kallisto/kallisto-0.43.0_1KGRef_Gencode19_k31/kallisto-0.43.0_1KGRef_Gencode19_k31.noGenes.index
GENOME_STAR_INDEX		$GENOME_STAR_INDEX_200

GENE_MODELS			${databaseDirectory}/gencode/gencode19/gencode.v19.annotation_plain.gtf
GENE_MODELS_EXCLUDE		${databaseDirectory}/gencode/gencode19/gencode.v19.annotation_plain.chrXYMT.rRNA.tRNA.gtf
GENE_MODELS_DEXSEQ		${databaseDirectory}/gencode/gencode19/gencode.v19.annotation_plain.dexseq.gff
GENE_MODELS_GC			${databaseDirectory}/gencode/gencode19/gencode.v19.annotation_plain.transcripts.autosomal_transcriptTypeProteinCoding_nonPseudo.1KGRef.gc
ARRIBA_KNOWN_FUSIONS		${hg19BaseDirectory}/tools_data/arriba/known_fusions_CancerGeneCensus_gencode19_2017-01-16.tsv.gz
ARRIBA_BLACKLIST		${hg19BaseDirectory}/tools_data/arriba/blacklist_hs37d5_gencode19_2017-01-09.tsv.gz
``
== Changelist

Release_1.0.11 [12th June 2017]

- Added README
- Added tmp-directory parameter to sambamba markdup (it previously wrote to /temp, but now it write to SCRATCH)
- Changed default scratch to outputAnalysisDir

Release_1.0.10 [9th June 2017]

- Added disableDoC_GATK flag for some samples where DoC calculations fail in RNAseQC
- Added additional flags for runnings kallisto in strand sepcific mode
- Corrected genome index for kallisto (perviously it used gene and transcript defintions from the GTF file, now it correctly only uses transcripts)

Release_1.0.8

- First official realease to OTP
