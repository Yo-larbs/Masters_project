#!/bin/bash
#please Specify whether virus or vector when running
if [ -z "$1" ]; then
  echo "Error: Please specify 'virus' or 'vector' when running the script."
  exit 1
fi

if [ "$1" != "virus" ] && [ "$1" != "vector" ]; then
  echo "Error: The argument must be either 'virus' or 'vector'."
  exit 1
fi
echo "Running PLINK2 VCF conversion for: $1"
Docker run --platform linux/amd64 -v ~/Documents/Masters_project:/data quay.io/biocontainers/plink2:2.00a2.3--hf22980b_0 plink2 --vcf "/data/output/VCF/${1}/Filtered_merged_genome/filtered_merged_genome.vcf.gz" --make-pgen --allow-extra-chr --out "/data/output/VCF/Plink_output/${1}/${1}_VCF"
Docker run --platform linux/amd64 -v ~/Documents/Masters_project:/data quay.io/biocontainers/plink2:2.00a2.3--hf22980b_0 plink2 --pfile  "/data/output/VCF/Plink_output/${1}/${1}_VCF" --pheno "/data/inputs/${1}_phenotype.txt" --glm firth --allow-extra-chr --pheno-name phenotype --out "/data/output/VCF/Plink_output/${1}/${1}_assoc_results"
