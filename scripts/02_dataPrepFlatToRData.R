# About -------------------------------------------------------------------
# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt
#
# 01/03/2019
# VCFs converted to flat files, now convert to RData to load for ShinyApp
# - input: flat files: "data\ShinyAppInput"
# - output: RData for shiny
# R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"

# Workspace ---------------------------------------------------------------
setwd("C:/Users/tdadaev/Desktop/Work/GitHubProjects/NGS")

library(data.table)

convertGT <- function(x){
  ifelse(x == "0/0", 0, 
         ifelse(x %in% c("0/1", "1/0"), 1,
                ifelse(x == "1/1", 2, 9)))
  # dots and half calls: set to missing:
  #    c(".", "./.", "0/.", "1/.") is  9
}

# Data --------------------------------------------------------------------
namesVCF <- gsub("_hg19.*", "", list.files("data/ShinyAppInput/", pattern = "*.GT"))


# ~Genotype ----------------------------------------------------------------
GT <- rbindlist(
  lapply(list.files("data/ShinyAppInput/", pattern = "*.GT", full.names = TRUE), function(i){
    #i = "data/ShinyAppInput/Familial_Exomes_hg19_VEPreannotated_filtered_sample_gencodev29cds10bpflank.vcf.gz.GT"
    s <- fread(gsub(".GT", ".SAMPLE", i, fixed = TRUE), header = FALSE)[, V1]
    d <- fread(i, drop = 1:5, col.names = s)
    cols <- colnames(d)
    d[ , (cols) := lapply(.SD, convertGT), .SDcols = cols ]
    
    varname <- fread(i, select = 1:5)
    varname <- paste(varname$V1, varname$V2, varname$V4, varname$V5, sep = "_" )
    
    melt(cbind(varname = varname, d),
         id.vars = "varname", variable.name = "sampleID")[ (value != 9), ]
    # res$vcf <- gsub("_hg19.*", "", basename(i))
    # res
  }))

# GT[ , cnt := .N, by = .(varname, sampleID) ]
# GT[ , cntUnique := length(unique(value)), by = .(varname, sampleID) ]
# xx <- GT[ cnt > 1, ]
# dim(xx)

# > dim(xx)
# [1] 391724      4
# > length(unique(xx$varname))
# [1] 69490
# > length(unique(xx$sampleID))
# [1] 774

# [1] 111108      4
# length(unique(xx$varname))
# [1] 1375
# length(unique(xx$sampleID))
# [1] 771
# GT$rn <- rowidv(GT, cols = c("varname", "sampleID"))
# GT[, .(rn := seq_len(.N)), by = .(varname, sampleID)]
# GT[ , rn := rowid(), by = .(varname, sampleID) ]
# table(GT[ !(G %in% c(0,1,2)), G])
# GT[ is.na(G), ]
# length(unique(GT$varname))

# to-do:get max when discordant, temp solution, will need to change...
GT[ , value := max(value), by = .(varname, sampleID) ]
GT <- unique(GT)
GT <- dcast(GT, varname ~ sampleID)

# fix MAF, flip: 0=common, 2=rare
GTalleleCounts <- data.table(
  AA = rowSums(GT[, -1] == 0, na.rm = TRUE),
  AB = rowSums(GT[, -1] == 1, na.rm = TRUE),
  BB = rowSums(GT[, -1] == 2, na.rm = TRUE),
  NN = rowSums(is.na(GT[, -1])))

# z <- x$AA + x$AB + x$BB
# xx <- data.table(
#   freqAA = x$AA/z,
#   freqAB = x$AB/z,
#   freqBB = x$BB/z)
# boxplot(xx)




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
  cbind(
    varname = paste(CLNSIG$CHROM, CLNSIG$POS, CLNSIG$REF, CLNSIG$ALT, sep = "_"),
    CLNSIG, CSQ)
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
           #i = list.files("data/GeneLists", full.names = TRUE, recursive = TRUE)[1]
           cbind(unique(fread(i, skip = 2, header = FALSE, col.names = "gene")),
                 panel = fread(i, nrows = 1, header = FALSE)[, V1], 
                 panelType = unlist(strsplit(i, "/", fixed = TRUE))[3])
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
pheno <- fread("data/20190315_progeny.csv",
               check.names = TRUE, na.strings = c("U", "", "NA"))
pheno <- pheno[ Study.ID %in% unique(unlist(SAMPLE)), ]

pheno[ , rs138213197 := factor(rs138213197) ]
pheno[ , EthnicityOA := factor(Onco_GenoAncestry) ]
pheno[ , PrCa := factor(Cancer.Confirmed.and.Unconfirmed.Prostate, levels = c(0, 1)) ]
pheno[ , COD_PrCa := factor(ifelse(tolower(Cause.of.death.is.PrCa) == "yes", 1,
                                   ifelse(tolower(Cause.of.death.is.PrCa) == "no", 0, NA))) ]
pheno[ , TStage := factor(gsub("T", "", TStage, fixed = TRUE)) ]
pheno[ , NStage := factor(gsub("N", "", NStage, fixed = TRUE)) ]
pheno[ , MStage := factor(gsub("M", "", MStage, fixed = TRUE)) ]

pheno[ , FH := factor(ifelse(FH == "Yes", 1, ifelse(FH == "No", 0, NA))) ]
pheno[ , PSADiag := round(PSADiag, 1) ]

pheno[ , NCCN := factor(NCCN, levels = c("Low", "Intermediate", "High", "VeryHigh", "Metastatic")) ]
pheno[ , NICE := factor(NICE, levels = c("Low", "Intermediate", "High")) ]

# add IDs that are not in Progeny (mostly controls?)
pheno <- merge(pheno, data.table(Study.ID = setdiff(colnames(GT[, -1]), pheno$Study.ID)), all = TRUE)


# Concordance -------------------------------------------------------------
GTc <- 
  lapply(list.files("data/ShinyAppInput/", pattern = "*.GT", full.names = TRUE), function(i){
    s <- fread(gsub(".GT", ".SAMPLE", i, fixed = TRUE), header = FALSE)[, V1]
    d <- fread(i, drop = 1:5, col.names = s)
    cols <- colnames(d)
    d[ , (cols) := lapply(.SD, convertGT), .SDcols = cols ]
    
    varname <- fread(i, select = 1:5)
    varname <- paste(varname$V1, varname$V2, varname$V4, varname$V5, sep = "_" )
    
    cbind(varname = varname, d)
  })
names(GTc) <- namesVCF

cc <- combn(namesVCF, m = 2)
overlapS <- apply(cc, 2, function(i){
  intersect(SAMPLE[[ i[1] ]], SAMPLE[[ i[2] ]])
  })
names(overlapS) <- apply(cc, 2, paste, collapse = "_and_")
# length(sort(unique(unlist(overlapS))))
# 1193

overlapV <- apply(cc, 2, function(i){
  intersect(ANNOT[[ i[1] ]]$varname, ANNOT[[ i[2] ]]$varname)
})
names(overlapV) <- apply(cc, 2, paste, collapse = "_and_")
# length(sort(unique(unlist(overlapV))))
# [1] 70461

# subset only sampels and vars that have any overlap with any other VCF
GTc <- lapply(GTc, function(i){
  # drop samples that have no overlapping vars even if samples overlap.
  samples <- unique(unlist(overlapS[ lengths(overlapV) > 0 ]))
  cols <- c("varname", intersect(colnames(i), samples))
  i[ varname %in% unique(unlist(overlapV)), ..cols]  
})
# lapply(GTc, dim)

# output ------------------------------------------------------------------
save(ANNOT, GT, GTalleleCounts, GTc, SAMPLE, geneListVCF, geneListPanel,
     namesVCF, filterCol, pheno,
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
