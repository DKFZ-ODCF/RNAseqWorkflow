## **Script to prepare gencode annotation data for RNA-seq analysis**

The reference data for GRCh38 and GRCm39 are prepared using the [refmake workflow](https://odcf-gitlab.dkfz.de/ODCF/refmake). This includes the
1. STAR index
2. Kallisto index
3. Gencode annotation GTF file

The downstream annotation files that are listed below are generated using the `prepare_gencode_annotation.sh` script. 
1. annotation.bed
2. annotation.nogene.gtf
3. annotation.chrXYMT.rRNA.gtf
4. annotation.dexseq.gff

The script can be run as follows:

```bash
sh prepare_gencode_annotation.sh /omics/odcf/reference_data/legacy/ngs_share/assemblies/hg_GRCh38/databases/gencode/GRCh38_decoy_ebv_alt_hla_phiX/gencode_v39_chr_patch_hapl_scaff/annotation.gtf
```

The Python script `dexseq_prepare_annotation2.py` was downloaded from [here](https://raw.githubusercontent.com/vivekbhr/Subread_to_DEXSeq/master/dexseq_prepare_annotation2.py) and edited to not print transcript IDs in the output files. This was to avoid the memory issue caused by the long concatenation of the transcript IDs.