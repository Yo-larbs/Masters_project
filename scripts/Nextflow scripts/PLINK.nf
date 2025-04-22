process PLINK {
    tag "PLINK"
    publishDir "${params.outdir}/Plink_output/${params.celltype}", mode:'copy'
    container "${ workflow.containerEngine == 'docker' && !task.ext.singularity_pull_docker_container ?
    'quay.io/biocontainers/plink2:2.00a2.3--hf22980b_0' :
    'biocontainers/plink2:2.00a2.3--hf22980b_0' }"
    input:
    path filteredmergedVCF
    path indexedmergedVCF
    path phenotype_file
    output:
    path "${params.celltype}_assoc_results.phenotype.glm.firth", emit: plink_results
    path "${params.celltype}_assoc_results.log", emit: plink_log
    script:
    """
    plink2 --vcf "$filteredmergedVCF" --make-pgen --allow-extra-chr --out "${params.celltype}_VCF"
    plink2 --pfile  "${params.celltype}_VCF" --pheno "${phenotype_file}" --glm firth --allow-extra-chr --pheno-name phenotype --out "${params.celltype}_assoc_results"
    """
}