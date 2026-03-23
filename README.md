# HIV masters project

## introduction

This is my Msc Dissertation project, i detailed my methods and code for posterity and for reproducibility and scrutiny. 

### experimental design

<img width="644" height="469" alt="image" src="https://github.com/user-attachments/assets/7a77ff48-8395-4975-bb0c-9040535fb003" />



HIV interaction with host genomic factors has been shown to drive HIV resistance in cells and possibly heterogeneity in infection. In this study 150 iPSCs were screened using a lentivirus (early events) and a HIV-1 strain (late events). These cell lines were then classified using PCA and clustering based on their response to infection. 

22 extremely susceptible and 25 extremely resistant cell lines classified in vector (early events) assay and 18 extremely susceptible and resistant cell lines classified in virus (late events) assay. Genomic association analysis was performed on all cells with PC1 from classification to measure phenotype using Plink2 and REGENIE, 

Differential expression analysis also performed using DEseq2 and RNA-seq of 54 of the 150 cell lines

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
vector/virus_phenotype.txt (categorical_only extremes)

vector/virus_plink_phenotype.txt (Quantitative using PC1)

### **NOTE** 
BILX_1 cell line has an typo (called BLIX_1) within the input file, i recommend manually changing within the phenotype file if you want to reporoduce


example Graphs produced:

<img width="438" height="313" alt="image" src="https://github.com/user-attachments/assets/5b17d486-ace2-40c4-a286-56f3868e787f" />


## HIV extreme comparison

### files 
HIV_vector.qmd

HIV_virus.qmd

### instructions
Need to have the vcfs for extreme phenotypes listed from cell line classification

annovar avinput file made from file is then used in annovar

**need to have annovar** (annovar code listed within the virus file)

### inputs: 
extreme VCFs (susceptible and resistant)

### outputs
avinput file (only comparing extremes)

example graph produced:

<img width="313" height="212" alt="image" src="https://github.com/user-attachments/assets/cd118bef-7dbe-4719-b1fb-828951d6fdf5" />


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

would also need  **annovar output** for all cell lines as well as **genotype file** (code used to produce included in Total_plink_virus_analysis)

example Graphs produced:



<img width="367" height="151" alt="image" src="https://github.com/user-attachments/assets/072c18c7-a682-499c-b8b6-0a7da4b1861d" />

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

example graph produced

<img width="348" height="261" alt="image" src="https://github.com/user-attachments/assets/aa54a7d0-5211-423a-a15c-ab2c543067d3" />



## eQTL

### files for results
eQTL done with PLink script
Plink_for_eQTL.sh

### files for analysis
Plink_eQTL.qmd is quarto documenting this

### NOTES
requires RNA counts (see transcriptomic), Annovar and genotype (see Plink)
exports TPM counts which is used as input in Plink_for_eQTL.sh


# SHINY APP

I made a shiny app to visualise the cleaning and classification sections of my work, just place the original data files in the same folder as the app (shiny-gene-list), Open the app using Rstudio and make sure you have Rshiny packages installed, if confused please consult comments in the app



