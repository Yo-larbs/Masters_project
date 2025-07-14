#!/bin/bash
#SBATCH --job-name=merge_batch
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=0-7    # Adjust based on how many batches you have
#SBATCH --cpus-per-task=8
export WORKDIR="/scratch/prj/ipsc/recovered/HIV_iPSC/Yoofi_project/Masters"
export EXOME_BED="$WORKDIR/GEN_exomic.bed"
 export Hg38="$WORKDIR/Hg38.fasta"
export IMG="$WORKDIR/bcftools:1.21.sif"
export SRC1="/scratch/prj/ipsc/recovered/HIV_iPSC/HipSci_Managed_Access/Exome/gVCF"
export SRC2="/scratch/prj/ipsc/recovered/HIV_iPSC/HipSci_Genome/gVCF"


find "$SRC1" "$SRC2" -name '*.g.vcf' > "$WORKDIR/all_gvcfs.txt"
split -l 20 all_gvcfs.txt gvcf_batch_

BATCH_FILES=($(ls gvcf_batch_*))  # Array of batch files
BATCH=${BATCH_FILES[$SLURM_ARRAY_TASK_ID]}

MERGED_BATCH_FILES=($(ls merged_gvcf_batch_*.vcf.gz))
MERGED_BATCH=${MERGED_BATCH_FILES[$SLURM_ARRAY_TASK_ID]}

singularity exec --bind "$WORKDIR":"$WORKDIR" "$IMG" bcftools merge -g $Hg38 -W --threads 8 -O z -o "$WORKDIR/merged_${BATCH##*/}.vcf.gz" $(cat "$BATCH")
singularity exec --bind "$WORKDIR":"$WORKDIR" "$IMG" bcftools view --threads 8 -W -R "$EXOME_BED" -O z -o "$WORKDIR/filtered_${MERGED_BATCH##*/}.vcf.gz" "$MERGED_BATCH"


FILTERED_MERGED_BATCH_FILES=($(ls filtered_merged_gvcf_batch_*.vcf.gz))
FILTERED_MERGED_BATCH=${FILTERED_MERGED_BATCH_FILES[$SLURM_ARRAY_TASK_ID]}

singularity exec --bind "$WORKDIR":"$WORKDIR" "$IMG" bcftools norm --threads 8 -W -f "$Hg38"  -O z -o "$WORKDIR/normalised_${MERGED_BATCH##*/}.vcf.gz"  "$FILTERED_MERGED_BATCH"

NORMALISED_MERGED_BATCH_FILES=($(ls normalised_merged_gvcf_batch_*.vcf.gz))
NORMALISED_MERGED_BATCH=${NORMALISED_MERGED_BATCH_FILES[$SLURM_ARRAY_TASK_ID]}

singularity exec --bind "$WORKDIR":"$WORKDIR" "$IMG" bcftools sort  -W -O z -o "$WORKDIR/sorted_${MERGED_BATCH##*/}.vcf.gz" "$NORMALISED_MERGED_BATCH"

SORTED_BATCH_FILES=($(ls sorted_merged_gvcf_batch_*.vcf.gz))  # Array of batch files
SORTED_BATCH=${SORTED_BATCH_FILES[$SLURM_ARRAY_TASK_ID]}
singularity exec --bind "$WORKDIR":"$WORKDIR" "$IMG" bcftools merge -g $Hg38 -W --threads 8 -O z -o "$WORKDIR/total_merged_${BATCH##*/}.vcf.gz" $(cat "$SORTED_BATCH")
