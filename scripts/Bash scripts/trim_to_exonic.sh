#!/bin/bash

Docker run --platform linux/amd64 -v ../..:/data quay.io/biocontainers/bcftools:1.7--0 bcftools view -R $1 $2 -o data/output/filtered_merged_genome.vcf && echo work_done
