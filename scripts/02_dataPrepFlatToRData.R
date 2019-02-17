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
    ifelse(x == "0/0", 0, 
           ifelse( x %in% c("0/1", "1/0"), 1,
                   ifelse(x == "1/1", 2, 9)))
    # dots and half calls: set to missing:
    #    c(".", "./.", "0/.", "1/.") is  9
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
  CLNSIG <- fread(i, sep = "\t", na.strings = c("", ".", "NA"))
  CLNSIG$CHROM <- as.character(CLNSIG$CHROM)
  #vep
  CSQ <- fread(gsub("CLNSIG", "CSQ", i, fixed = TRUE), sep = "|", na.strings = c("", ".", "NA"))
  cbind(CLNSIG, CSQ)
})
names(ANNOT) <- namesVCF


# ~GeneSymbols VCF --------------------------------------------------------
#geneListVCF <- sort(unique(unlist(lapply(ANNOT, function(i) i[ , SYMBOL]))))
geneListVCF <- rbindlist(
  lapply(names(ANNOT), function(i) {
    data.table(gene = sort(unique(ANNOT[[ i ]][ SYMBOL != "", SYMBOL])),
               panel = i)
  }))

# ~GeneSymbols Panel ------------------------------------------------------
geneListPanel <- rbindlist(
  lapply(list.files("data/GeneLists", full.names = TRUE, recursive = TRUE),
         function(i){
           #i = list.files("data/GeneLists/", full.names = TRUE, recursive = TRUE)[1]
           cbind(unique(fread(i, skip = 2, header = FALSE, col.names = "gene")),
                 panel = fread(i, nrows = 1, header = FALSE)[, V1])
           
         }))

# ~Annot filters -----------------------------------------------------------
x <- unique(rbindlist(ANNOT)[, c("CLNSIG", "Consequence", "IMPACT",
                                 "SIFT", "PolyPhen", "LoF")])
filterCol <- list(
  CLNSIG = sort(unique(gsub("^_", "", unlist(strsplit(x$CLNSIG[ !is.na(x$CLNSIG) ], "[,/]"))))),
  Consequence = sort(unique(unlist(strsplit(x$Consequence, "&")))),
  IMPACT = sort(unique(x$IMPACT[ !is.na(x$IMPACT) ])),
  SIFT = sort(unique(gsub("\\(.*", "", x$SIFT[ !is.na(x$SIFT) ]))),
  PolyPhen = sort(unique(gsub("\\(.*", "", x$PolyPhen[ !is.na(x$PolyPhen) ]))),
  LoF = sort(unique(gsub("\\(.*", "", x$LoF[ !is.na(x$LoF) ]))))



# ~ Phenotype ---------------------------------------------------------------
pheno <- fread("data/20190215_progeny.csv",
               check.names = TRUE, na.strings = c("U", "", "NA"))
pheno <- pheno[ Study.ID %in% unique(unlist(SAMPLE)), ]

# output ------------------------------------------------------------------
save(ANNOT, GT, SAMPLE, geneListVCF, geneListPanel, namesVCF, filterCol, pheno,
     file = "data/data.RData")


# Testing -----------------------------------------------------------------
# GTsummary <- lapply(list.files("data/ShinyAppInput/", pattern = "*.GT", full.names = TRUE), function(i){
#   d <- fread(i)
#   d <- d[ , 6:ncol(d)]
#   x <- table(unlist(d))
#   res <- merge(data.frame(x), data.frame(round(prop.table(x), 2)), by = "Var1")
#   colnames(res) <- c("GT", "N", "%")
#   res
# })
# names(GTsummary) <- list.files("data/ShinyAppInput/", pattern = "*.GT")
