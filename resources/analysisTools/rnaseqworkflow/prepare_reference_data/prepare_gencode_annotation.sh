#! /bin/env bash
# Prepare gencode annotation for RNAseq workflow based on the GTF file from the refmake workflow
# Example run sh prepare_gencode_annotation.sh /omics/odcf/reference_data/legacy/ngs_share/assemblies/hg_GRCh38/databases/gencode/GRCh38_decoy_ebv_alt_hla_phiX/gencode_v39_chr_patch_hapl_scaff/annotation.gtf
module load bedops/2.4.14
GTF_FILE=$1
GTF_FILE_basename=${GTF_FILE%.gtf}

# GTF to bed
convert2bed -i gtf < $GTF_FILE > ${GTF_FILE_basename}.bed

## For RNAseq workflow
cat $GTF_FILE | grep -P "^chr[X|Y|MT]|\brRNA" | awk -F '\t' '$3=="gene"' > ${GTF_FILE_basename}.chrXYMT.rRNA.gtf

# Creating dexseq gff files using vivekbhr/Subread_to_DEXSeq
# wget https://raw.githubusercontent.com/vivekbhr/Subread_to_DEXSeq/master/dexseq_prepare_annotation2.py
# conda create --name htseq -c bioconda htseq
source activate htseq

python dexseq_prepare_annotation2.py -f ${GTF_FILE_basename}.dexseq.gtf $GTF_FILE ${GTF_FILE_basename}.dexseq.gff

## Creating nogene GTF file for RNAseqQC analysis
# The gene annotations doesn't contain 'transcript_id', so has to be excluded for the RNAseqQC analysis
cat $GTF_FILE | awk '$3 != "gene"' > ${GTF_FILE_basename}.nogene.gtf
