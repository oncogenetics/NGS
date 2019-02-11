# About -------------------------------------------------------------------
# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt
#
# 08/02/2019
# VCFs converted to flat files, now convert to RData to load for ShinyApp
# - input: flat files: "data\ShinyAppInput"
# - output: RData for shiny
# R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"

# Workspace ---------------------------------------------------------------
setwd("C:/Users/tdadaev/Desktop/Work/GitHubProjects/NGS")

library(data.table)


# Data --------------------------------------------------------------------
namesVCF <- gsub("_hg19.*", "", list.files("data/ShinyAppInput/", pattern = "*.GT"))


# ~Genotype ----------------------------------------------------------------
GT <- lapply(list.files("data/ShinyAppInput/", pattern = "*.GT", full.names = TRUE), function(i){
  #i = "data/ShinyAppInput/Familial_Exomes_hg19_VEPreannotated_filtered_sample_gencodev29cds10bpflank.vcf.gz.GT"
  d <- fread(i)
  d <- d[ , 6:ncol(d)]
  d[] <- lapply(d, function(x){
    ifelse(x %in% c("0","0 0","0/0"), 0, 
           ifelse( x %in% c("0/1", "1 0", "1/0"), 1,
                   ifelse(x == "1/1", 2, 9))) 
    # . dots and half calls: set to missing
    # c(".", "./.", "0/.", "1/.") is  9
  })
  d <- as.matrix(d)
  colnames(d) <- NULL
  d
})
names(GT) <- namesVCF

# ~Sample ID ---------------------------------------------------------------
SAMPLE <- lapply(list.files("data/ShinyAppInput/", pattern = "*.SAMPLE", full.names = TRUE), function(i){
  fread(i, header = FALSE)[, V1]
  })
names(SAMPLE) <- namesVCF



# ~Annot: clinvar+vep ------------------------------------------------------
ANNOT <- lapply(list.files("data/ShinyAppInput/", pattern = "*.CLNSIG", full.names = TRUE), function(i){
  #i = "data/ShinyAppInput/Familial_Exomes_hg19_VEPreannotated_filtered_sample_gencodev29cds10bpflank.vcf.gz.CLNSIG"
  #clinvar
  CLNSIG <- fread(i)
  #vep
  CSQ <- fread(gsub("CLNSIG", "CSQ", i, fixed = TRUE))
  cbind(CLNSIG, CSQ)
})
names(ANNOT) <- namesVCF

# output ------------------------------------------------------------------
save(ANNOT, GT, SAMPLE, namesVCF, file = "data/data.RData")




# Testing -----------------------------------------------------------------
# GTsummary <- lapply(list.files("data/ShinyAppInput/", pattern = "*.GT", full.names = TRUE), function(i){
#     #i = "data/ShinyAppInput/Familial_Exomes_hg19_VEPreannotated_filtered_sample_gencodev29cds10bpflank.vcf.gz.GT"
#   d <- fread(i)
#   d <- d[ , 6:ncol(d)]
#   x <- table(unlist(d))
#   res <- merge(data.frame(x), data.frame(round(prop.table(x), 2)), by = "Var1")
#   colnames(res) <- c("GT", "N", "%")
#   res
# })
# names(GTsummary) <- list.files("data/ShinyAppInput/", pattern = "*.GT")
# $AEPv2_hg19_VEPreannotated_filtered_samplerenamed_gencodev29cds10bpflank.vcf.gz.GT
# GT        N    %
# 1 ./.  1261765 0.02
# 2 0/0 44752969 0.85
# 3 0/1  4080592 0.08
# 4 1/1  2541314 0.05
# 
# $DRG_2441_hg19_VEPreannotated_filtered_samplerenamed_gencodev29cds10bpflank.vcf.gz.GT
# GT        N    %
# 1   .   124815 0.01
# 2 0/.      190 0.00
# 3 0/0 13673109 0.97
# 4 0/1   204600 0.01
# 5 1/.        1 0.00
# 6 1/1    94060 0.01
# 
# $Eeles_BRCA1_Sanger_hg19_VEPreannotated_150119.vcf.gz.GT
# GT    N %
# 1   0    4 0
# 2 0 0 3520 1
# 3 1 0    4 0
# 
# $Eeles_BRCA2_FCAP_hg19_VEPreannotated_150119.vcf.gz.GT
# GT     N %
# 1   0    20 0
# 2 0 0 36518 1
# 3 1 0    22 0
# 
# $Familial_Exomes_hg19_VEPreannotated_filtered_sample_gencodev29cds10bpflank.vcf.gz.GT
# GT       N    %
# 1 ./.   26687 0.01
# 2 0/0 3679495 0.74
# 3 0/1  738366 0.15
# 4 1/0    1967 0.00
# 5 1/1  493867 0.10
# 
# $MCK_191_hg19_VEPreannotated_filtered_sample_gencodev29cds10bpflank.vcf.gz.GT
# GT     N    %
# 1 ./.   315 0.00
# 2 0/0 67400 0.90
# 3 0/1  4567 0.06
# 4 1/0    53 0.00
# 5 1/1  2346 0.03




# vcfR didn't work, too heavy on memory
# library(vcfR)
# 
# fileVCF_Familial ="data/DRG_2441_hg19_VEPreannotated_filtered_samplerenamed_gencodev29cds10bpflank.vcf.gz"
# vcf <- read.vcfR(fileVCF_Familial)#, verbose = FALSE )
# 
# chromoqc(vcf, dp.alpha=20)
# 
# vcf@fix
# x <- vcfR2tidy(vcf)
