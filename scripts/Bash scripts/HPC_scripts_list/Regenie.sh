#!/bin/bash
export WORKDIR="/scratch/prj/ipsc/recovered/HIV_iPSC/Yoofi_project/Masters"
export PLINK2="$WORKDIR/plink2.sif"
export REGENIE="$WORKDIR/regenie.sif"

#awk 'BEGIN {OFS="\t"; print "\#FID","IID","SEX"} NR>1 {OFS="\t"; print $1, $1, "0"}' "$WORKDIR/pheno_plink.txt" > tmp.txt


#Regenie uses plink2 output file here i generate the PSAM (sample file), Pvar and Pgen (genotype)

singularity exec --bind "$WORKDIR":"$WORKDIR" "$PLINK2" plink2 --vcf "$WORKDIR/all_batches_merged/total_merged.normalized.vcf.gz"  --make-pgen --chr chr1-chr22 --mac 5 --maf 0.0001  --set-all-var-ids @:#:\$r:\$a  --rm-dup force-first  --allow-extra-chr --new-id-max-allele-len 393 --out "$WORKDIR/plink_files/total_plink_VCF"

#Psam is only 2 column output, i need it to have 3 columns FID (family ID) IID (individual ID) and sex
#becuase using individual iPSC, FID and IID would be the same

awk 'BEGIN {OFS="\t"; print "\#FID","IID","SEX"} NR>1 {OFS="\t"; print $1, $1, "0"}' "$WORKDIR/plink_files/total_plink_VCF.psam" > tmp.txt
rm "$WORKDIR/plink_files/total_plink_VCF.psam"
mv tmp.txt "$WORKDIR/plink_files/total_plink_VCF.psam"

# here we get the frequency of all variants (this was previously to filter out but not needed anymore 
#singularity exec --bind "$WORKDIR":"$WORKDIR" "$PLINK2" plink2 --pfile "$WORKDIR/plink_files/total_plink_VCF" --freq \
#  --out "$WORKDIR/plink_files/total_plink_VCF"

#awk '$5 == 0 { print $2 }' "$WORKDIR/plink_files/total_plink_VCF.afreq" > zero_variants.txt

# LD pruning 

singularity exec --bind "$WORKDIR":"$WORKDIR" "$PLINK2" plink2 \\
  --pfile "$WORKDIR/plink_files/total_plink_VCF" \
  --indep-pairwise 200 50 0.1 \
  --out "$WORKDIR/plink_files/total_plink_VCF"

#PCA

singularity exec --bind "$WORKDIR":"$WORKDIR" "$PLINK2" plink2 \
  --pfile "$WORKDIR/plink_files/total_plink_VCF" \
  --extract "$WORKDIR/plink_files/total_plink_VCF.prune.in" \
  --pca 4 \
  --out "$WORKDIR/plink_files/total_plink_VCF"

# Need to combine  the PCA into a covariate file we merge by using the IID

( printf "FID\tIID\tPHENO\tPC1\tPC2\tPC3\tPC4\n"
  join -t $'\t' -1 2 -2 2 \
       <( sort -k2,2 "$WORKDIR/pheno_vector_plink.txt" ) \
       <(tail -n +2 "$WORKDIR/plink_files/total_plink_VCF.eigenvec" | sort -k2,2  ) \
  | awk 'BEGIN{OFS="\t"} { print $2,$1,$3,$5,$6,$7,$8 }'
) > re_pheno_covar.txt

#Regenie step 1 feeding it pruned LD  in prediction
singularity exec --bind "$WORKDIR":"$WORKDIR" "$REGENIE" regenie \
  --step 1 \
  --pgen "$WORKDIR/plink_files/total_plink_VCF" \
  --bsize 1000 \
 --exclude "$WORKDIR/plink_files/total_plink_VCF.prune.out" \
  --threads 8 \
 --phenoFile pheno_vector_plink.txt  \
  --out "$WORKDIR/regenie_files/re_step1"

singularity exec --bind "$WORKDIR":"$WORKDIR" "$REGENIE" regenie \
 --step 2 \
 --bsize 1000 \
  --pgen "$WORKDIR/plink_files/total_plink_VCF" \
  --phenoFile re_pheno_covar.txt \
 --phenoCol PHENO \
 --minMAC 5 \
  --covarFile re_pheno_covar.txt --covarCol PC1,PC2,PC3,PC4 \
  --pred "$WORKDIR/regenie_files/re_step1_pred.list" \
  --threads 8 \
  --out "$WORKDIR/regenie_files/re_assoc"
