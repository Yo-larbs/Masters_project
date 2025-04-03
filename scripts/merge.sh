#!/bin/bash

Docker run --platform linux/amd64 -v /Users/yoofilarbi/Documents/Masters_project:/data  quay.io/biocontainers/bcftools:1.7--0 bcftools merge "data/output/reheadered/$(basename $1)" "data/output/reheadered/$(basename $2)" "data/output/reheadered/$(basename $3)"  -o /data/output/merged.vcf
Docker run --platform linux/amd64 -v /Users/yoofilarbi/Documents/Masters_project:/data quay.io/biocontainers/bcftools:1.7--0 bcftools convert -O z -o  "/data/output/merged.vcf.gz" "/data/output/merged.vcf" 
Docker run --platform linux/amd64 -v /Users/yoofilarbi/Documents/Masters_project:/data quay.io/biocontainers/bcftools:1.7--0 bcftools index  "/data/output/merged.vcf.gz" && rm "../output/merged.vcf"
