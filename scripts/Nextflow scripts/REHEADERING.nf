process REHEADERING {
    tag "${vcf_file.simpleName}"
    publishDir "${params.outdir}/${params.celltype}/reheadered", mode: 'copy'
    container "${ workflow.containerEngine == 'docker' && !task.ext.singularity_pull_docker_container ?
    'quay.io/biocontainers/bcftools:1.9--ha228f0b_4':
    'bcftools:1.9--ha228f0b_4' }"

    input:
        tuple val(fullpath), path(vcf_file)
    output:
        path "reheadered_*.vcf.gz", emit: reheaderedVCFs
        path "reheadered_*.vcf.gz.csi", emit: indexes
    script:

    """
    echo "$fullpath"
    SampleName=\$(basename "${vcf_file}" | sed -E 's/(\\.filtered\\..*\$|_ERR.*\$)//')
    echo "\$SampleName" > "\${SampleName}.txt"
    
    bcftools reheader -s "\${SampleName}.txt" "./\$(basename ${vcf_file})" -o "reheadered_\${SampleName}.vcf" && echo "work_done"
    rm "\${SampleName}.txt"
    
    bcftools convert -O z -o "reheadered_\${SampleName}.vcf.gz" "reheadered_\${SampleName}.vcf"
    
    bcftools index "reheadered_\${SampleName}.vcf.gz" && rm "reheadered_\${SampleName}.vcf"
    """
}
