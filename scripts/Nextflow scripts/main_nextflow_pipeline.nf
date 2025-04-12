#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
params.vcfPattern = "${params.projectDir}/inputs/VCF/Vector/**/*.vcf"
params.bedfile= "${params.projectDir}/inputs/exomic.bed"

include {REHEADERING} from './REHEADERING'

include {MERGE} from './MERGED'

include {TRIMMING} from './TRIMMING'




workflow{
    Channel.fromPath(params.vcfPattern)
    .map { file -> tuple(file, file.toAbsolutePath().toString()) }
    .set{vcf_files}
    Channel.fromPath(params.bedfile).set{bedfile}
    REHEADERING(vcf_files)
    MERGE(REHEADERING.out.reheaderedVCFs.collect(),REHEADERING.out.indexes.collect())
    TRIMMING(bedfile, MERGE.out)
}