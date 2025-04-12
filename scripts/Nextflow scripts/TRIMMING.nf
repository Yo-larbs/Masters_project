process TRIMMING{
    tag "${mergedVCF.simpleName}"
    publishDir "${params.outdir}/Filtered_merged_genome", mode: 'copy'
    container "${ workflow.containerEngine == 'docker' && !task.ext.singularity_pull_docker_container ?
    'quay.io/biocontainers/bcftools:1.7--0' :
    'biocontainers/bcftools:1.7--0' }"
    input:
    path bedfile 
    path mergedVCF
    path indexedmergedVCF
    output:
    path "filtered_merged_genome.vcf", emit: filtered_merged
    path "filtered_merged_genome.vcf.gz", emit: filteredmergedVCF
    path "filtered_merged_genome.vcf.gz.csi", emit: indexedmergedVCF

    script:
    """
    bcftools view -R $bedfile $mergedVCF -o "filtered_merged_genome.vcf" && echo work_done
    bcftools sort -o "filtered_merged_genome.vcf" "filtered_merged_genome.vcf"
    bcftools convert -O z -o  "filtered_merged_genome.vcf.gz" "filtered_merged_genome.vcf"
    bcftools index  "filtered_merged_genome.vcf.gz" 
    """
}