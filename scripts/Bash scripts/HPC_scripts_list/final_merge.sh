#!/bin/bash

export WORKDIR="/scratch/prj/ipsc/recovered/HIV_iPSC/Yoofi_project/Masters"
export EXOME_BED="$WORKDIR/GEN_exomic.bed"
 export Hg38="$WORKDIR/Hg38.fasta"
export IMG="$WORKDIR/bcftools:1.21.sif"

#merge all batch files

SORTED_BATCH_FILES=$(ls sorted_merged_gvcf_batch_*.vcf.gz)  # Array of batch files
singularity exec --bind "$WORKDIR":"$WORKDIR" "$IMG" bcftools merge -g $Hg38 -W --threads 12 -O z -o "$WORKDIR/total_merged_gvcf.vcf.gz" $SORTED_BATCH_FILES

# get rid off large structural and complex changes and any spot that has no mutation
singularity exec --bind "$WORKDIR":"$WORKDIR" "$IMG" bcftools view --threads 12 -W -v snps,indels -O z  -o "$WORKDIR/filtered_total_merged_gvcf.vcf.gz" "$WORKDIR/total_merged_gvcf.vcf.gz"

# get rid of multi allelics
singularity exec --bind "$WORKDIR":"$WORKDIR" "$IMG" bcftools norm -W -m -both -O z -o "$WORKDIR/total_merged.vcf.gz" "$WORKDIR/filtered_total_merged_gvcf.vcf.gz"

#get rid of NON_REF columns
singularity exec --bind "$WORKDIR":"$WORKDIR" "$IMG" bcftools view \
  -v snps,indels -W\
  -e 'ALT="<NON_REF>"' \
  -O z -o "$WORKDIR/total_merged.normalized.vcf.gz" \
  "$WORKDIR/total_merged.vcf.gz"
