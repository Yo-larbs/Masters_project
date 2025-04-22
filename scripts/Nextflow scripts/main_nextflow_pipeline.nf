#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
params.celltype= "vector"
params.vcfPattern = "${params.projectDir}/inputs/VCF/${params.celltype}/**/*.vcf"
params.bedfile= "${params.projectDir}/inputs/exomic.bed"
params.phenotype="${params.projectDir}/inputs/${params.celltype}_phenotype.txt"

include {REHEADERING} from './REHEADERING'

include {MERGE} from './MERGED'

include {TRIMMING} from './TRIMMING'

include {PLINK} from './PLINK'




workflow{
    Channel.fromPath(params.vcfPattern)
    .map { file -> tuple(file, file.toAbsolutePath().toString()) }
    .set{vcf_files}
    Channel.fromPath(params.bedfile).set{bedfile}
    Channel.fromPath(params.phenotype).set{phenotype_file}
    REHEADERING(vcf_files)
    MERGE(REHEADERING.out.reheaderedVCFs.collect(),REHEADERING.out.indexes.collect())
    TRIMMING(bedfile, MERGE.out)
    PLINK(TRIMMING.out.filteredmergedVCF,TRIMMING.out.indexedmergedVCF,phenotype_file)
}