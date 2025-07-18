process MERGE {
    tag "${reheaderedVCF.simpleName}"
    publishDir "${params.outdir}/${params.celltype}/Merged_genome", mode: 'copy'
    container "${ workflow.containerEngine == 'docker' && !task.ext.singularity_pull_docker_container ?
    'quay.io/biocontainers/bcftools:1.9--ha228f0b_4' :
    'biocontainers/bcftools:1.9--ha228f0b_4' }"

    input:
        path reheaderedVCF
        path indexes
    output:
        path "merged.vcf.gz", emit: mergedVCF
        path "merged.vcf.gz.csi", emit: indexedmergedVCF
    script:
    """
    bcftools merge $reheaderedVCF  -o "merged.vcf"
    bcftools convert -O z -o  "merged.vcf.gz" "merged.vcf"
    bcftools index  "merged.vcf.gz" && rm "merged.vcf"
    """
}
