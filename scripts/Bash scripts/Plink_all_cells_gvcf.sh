#!/bin/bash
#when running make sure VCF merged file and the phenotype file is in same folder as the script
#please use virus or vector as an arguments (no caps)
Docker run --rm --platform linux/amd64 -v "$PWD":/data quay.io/biocontainers/plink2:2.00a2.3--hf22980b_0 plink2 --vcf "/data/total_merged.normalized.vcf.gz" --make-pgen --allow-extra-chr  --out "/data/total_plink_VCF"
Docker run --rm --platform linux/amd64 -v "$PWD":/data quay.io/biocontainers/plink2:2.00a2.3--hf22980b_0 plink2 --pfile  "/data/total_plink_VCF" --mac 5 --pheno "/data/plink_${1}_phenotype.txt" --glm firth-fallback --allow-extra-chr --pheno-name PC1 --out "/data/total_${1}_assoc_results"
