#!/bin/bash
Docker run --platform linux/amd64 -v ../..:/data quay.io/biocontainers/bcftools:1.7--0 bcftools sort -o "/data/output/filtered_merged_genome.vcf" "/data/output/filtered_merged_genome.vcf"
Docker run --platform linux/amd64 -v ../..:/data quay.io/biocontainers/bcftools:1.7--0 bcftools convert -O z -o  "/data/output/filtered_merged_genome.vcf.gz" "/data/output/filtered_merged_genome.vcf"
Docker run --platform linux/amd64 -v ../..:/data quay.io/biocontainers/bcftools:1.7--0 bcftools index  "/data/output/filtered_merged_genome.vcf.gz"
