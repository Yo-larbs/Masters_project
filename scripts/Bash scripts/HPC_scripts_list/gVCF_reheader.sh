#!/bin/bash

module load parallel
export SRC1="/scratch/prj/ipsc/recovered/HIV_iPSC/HipSci_Managed_Access/Exome/gVCF"
export SRC2="/scratch/prj/ipsc/recovered/HIV_iPSC/HipSci_Genome/gVCF"
export WORKDIR="/scratch/prj/ipsc/recovered/HIV_iPSC/Yoofi_project/Masters"
export EXOME_BED="masters/GEN_exomic.bed"
 export Hg38="masters/hg38.fa"


g="$1"
SampleName=$(basename "$g" | sed -E 's/(\.g\..*$|_ERR.*$)//')

echo "$SampleName" > "$WORKDIR/sorted/${SampleName}.txt"

bcftools reheader --threads 12 -s "$WORKDIR/sorted/${SampleName}.txt" "$g" | bcftools view --threads 12 -O z -o "$WORKDIR/sorted/${SampleName}.reheadered.vcf.gz"
bcftools index "$WORKDIR/sorted/${SampleName}.reheadered.vcf.gz"

rm "$WORKDIR/sorted/${SampleName}.txt"

