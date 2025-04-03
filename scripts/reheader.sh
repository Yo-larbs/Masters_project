
#!/bin/bash

SampleName=$(basename "$1"|sed -E 's/(\.filtered\..*$|_ERR.*$)//')
#SampleName=${SampleName%_}
DATA="/Users/yoofilarbi/Documents/Masters_project/inputs/VCF/"
echo "$SampleName" > "${SampleName}.txt"

Docker run --platform linux/amd64 -v /Users/yoofilarbi/Documents/Masters_project:/data  quay.io/biocontainers/bcftools:1.7--0 bcftools reheader -s "/data/scripts/${SampleName}.txt" "/data/inputs/VCF/$(basename $1)" -o "/data/output/reheadered/reheadered_${SampleName}.vcf" && echo "work_done"
rm "${SampleName}.txt"

Docker run --platform linux/amd64 -v /Users/yoofilarbi/Documents/Masters_project:/data quay.io/biocontainers/bcftools:1.7--0 bcftools convert -O z -o "/data/output/reheadered/reheadered_${SampleName}.vcf.gz" "/data/output/reheadered/reheadered_${SampleName}.vcf"

Docker run --platform linux/amd64 -v /Users/yoofilarbi/Documents/Masters_project:/data quay.io/biocontainers/bcftools:1.7--0 bcftools index  "/data/output/reheadered/reheadered_${SampleName}.vcf.gz" && rm "../output/reheadered/reheadered_${SampleName}.vcf"
