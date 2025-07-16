# HIV masters project

## General workflow

1. cleaning and classification
2. genomic
3. plink/REGENIE
4. transcriptomic
5. eQTL

## Cleaning and classification
Vector (early events assay) - Cleaning_and_classification vector.rmd

Virus (later events assay) - Cleaning_and_classifcation virus.rmd

### inputs
hiv1_vectors_collated_n0_uncleaned(Sheet1).csv

hiv1_virus_colllated_n0_uncleaned(Sheet1).csv

### outputs
vector/virus_phenotype (categorical_only extremes)

vector/virus_plink_phenotype (Quantitative using PC1)

### **NOTE** 
BILX_1 cell line has an typo (called BLIX_1) within the input file, i recommend manually changing within the phenotype file if you want to reporoduce

## HIV extreme comparison

### files 
HIV_vector.qmd

HIV_virus.qmd

### instructions
Need to have the vcf for extremes listed from cell line classification
annovar avinput file made from file is then used in annovar
**need to have annovar** (annovar code listed within the virus file)

### inputs: 
extreme VCFs (susceptible and resistant)

## Plink2
combining all gVCF was done on Kings create HPC https://docs.er.kcl.ac.uk/CREATE/access/ .would need singularity pull files on to HPC, also uses specific file Paths and os would need some editing if trying to replicate

### files for gVCF merge
codes ran in order: gVCF_reheader.sh, gVCF_merge, final_merge.sh
egressed combined.VCF onto personal device

### files for results
Plink_all_cells_gvcf.sh
use either "vector" or "virus" as argument make sure Docker is on
should give .assoc file

### files for analysis
Total_plink_vector_analysis.qmd
Total_plink_virus_analysis.qmd

would also need  annovar output for all cell lines as well as genotype file (code used to produce included in Total_plink_virus_analysis)

## REGENIE

### files for results
run Regenie.sh , ran on HPC, using combined vcf

### files for analysis
Regenie_vector.qmd

Regenie_virus.qmd

### NOTE
these also need annovar and genotype file, but same ones as the ones documented in plink qmd can be used 

## Transcriptomic
Used Counts file, path specified in file,
contained 54/150 cell lines, multiple replicates for each

### files for analysis
RNA-seq_virus,
RNA-seq_vector


## eQTL

### files for results
eQTL done with PLink script
Plink_for_eQTL.sh

### files for analysis
Plink_eQTL.qmd is quarto documenting this

### NOTES
requires RNA counts (see transcriptomic), Annovar and genotype (see Plink)
exports TPM counts which is used as input in Plink_for_eQTL.sh




