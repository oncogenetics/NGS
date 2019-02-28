# About -------------------------------------------------------------------
# Adapted from Burcu's code at:
#/auto/pmd-02/wesp/exampleCode/burdenTesting.R
# - read GT per chrom, and run burden test
#   - Tier1
#   - Tier2


# Input -------------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

chrN = args[ 1 ]
chr = ifelse(chrN == 23, "X", ifelse(chrN == 24, "Y", as.character(chrN)))
chrom = paste0("chr", chr)
fileTier = args[ 2 ]
folderOutput = args[ 3 ]

# chrN=23
# chr="Y"
# chrom="chrY"
# fileTier="/auto/pmd-02/wesp/td_844/siteFiles/20181112/v_02_T1_repeat.txt"
# folderOutput="/auto/pmd-02/wesp/td_844/output/20181121_simpleBurden"
fileOutput = paste0(folderOutput, "/",
                    tools::file_path_sans_ext(basename(fileTier)), "_",
                    chrN,".txt")

# Workspace ---------------------------------------------------------------
.libPaths(c(.libPaths(), "/auto/pmd-02/wesp/software/R_3.5.0_libs"))

#setwd("C:/Users/tdadaev/Desktop/NIH/burden_Rglm")
setwd("/auto/pmd-02/wesp/td_844/burden_Rglm")

library("data.table")

#fileGenes = "../data/refFlat_hg19/refFlat_hg19_matchingGenes.txt"
fileGenes = "/auto/pmd-02/wesp/software/data/refFlat_hg19_chunks/refFlat_hg19_matchingGenes.txt"

fileGT = paste0("/auto/pmd-02/wesp/data/withGenotypeRefinement_chrom/", chrN, "_GT.txt")
fileAF = paste0("/auto/pmd-02/wesp/data/withGenotypeRefinement_chrom/", chrN, "_AF.txt")

fileCovar = "/auto/pmd-02/wesp/data/covariates/wesp_covar.covar"
filePheno = "/auto/pmd-02/wesp/data/covariates/wesp_pheno.pheno"  
fileSampleID = "/auto/pmd-02/wesp/data/covariates/wesp_id_VCF.txt"


# Data --------------------------------------------------------------------
# read in covariate & phenotype files
covar <- fread(fileCovar)
pheno <- fread(filePheno)
# all(covar$fid == pheno$iid)
# true
# all(covar$STATUS == pheno$STATUS)
# false

gt <- fread(fileGT)
af <- fread(fileAF, col.names = c("chr", "pos", "ref", "alt", "filter", "af"))
af$rowID <- seq(nrow(af))

# sampleID
setnames(gt, fread(fileSampleID, header = FALSE)[, V1])

# Tier vars: chr pos
tierChrPos <- fread(fileTier, col.names = c("chr", "pos"))
x = chr
tierChrPos <- tierChrPos[ chr == x, ]

genes <- fread(fileGenes,
               col.names = c("gene", "region", "chr", "start", "end", "width", "chunk"))
genes <- genes[ chr == chrom, ]

# match vars with gene position, find overlap
x1 <- af[ pos %in% tierChrPos$pos, list(start = pos, end = pos, rowID)]
x2 <- genes[, list(start, end, gene)]
setkeyv(x1, c("start", "end"))
setkeyv(x2, c("start", "end"))
tierGene <- foverlaps(x1, x2)
tierGene <- tierGene[ !is.na(gene), ]

# Burden Testing -----------------------------------------------------------
# take Status from pheno, make binary: 1,2 to 0,1
glmDat <- covar[, list(caco = pheno$STATUS - 1, age, PC1, PC2, PC3, STUDY)]

null.model = glm(caco ~ age + PC1 + PC2 + PC3 + STUDY, family = binomial, data = glmDat)
#summary(nullModel)

output <- rbindlist(
  lapply(unique(tierGene$gene), function(i){
    # i = "SYCE1"
    variantsRowID <- unique(tierGene[ gene == i, rowID])
    #if(length(variantsRowID) > 0) {
    G <- gt[ variantsRowID, covar$fid, with = FALSE ]
    
    #flip G based on MAF
    G[ (af[ variantsRowID, af] > 0.5), names(G) := lapply(.SD, function(x) 2L - x) ]

    # NA(9) or missing to 0(REF)
    G[, names(G) := lapply(.SD, function(x) ifelse(x %in% c(0:2), x, 0L)) ]

    glmDat[, X:=colSums(G) ]
    Xco = sum(glmDat[ caco == 0, X])
    Xca = sum(glmDat[ caco == 1, X])
    
    if(sum(glmDat$X) > 5){
      reg <- glm(caco ~ X + age + PC1 + PC2 + PC3 + STUDY, family = binomial, data = glmDat)
      s.reg <- summary(reg)
      
      chi.LRT <- 2 * (logLik(reg) - logLik(null.model))
      if(is.na(chi.LRT)){ P.LRT <- NA } else {
        P.LRT <- pchisq(chi.LRT, df = 1, lower.tail = FALSE) }
      
      Xest <- s.reg$coef["X", "Estimate"]
      Xse <- s.reg$coef["X", "Std. Error"]
      
      res <- data.frame(gene = i,
                        numVar = nrow(G),
                        allelesCo = Xco,
                        allelesCa = Xca,
                        OR = exp(Xest),
                        Lower.CI = exp(Xest - 1.96 * Xse),
                        Upper.CI = exp(Xest + 1.96 * Xse),
                        P = s.reg$coef["X", "Pr(>|z|)"],
                        chi.LRT = chi.LRT,
                        P.LRT = P.LRT)
    } else {
      res <- data.frame(gene = i,
                        numVar = nrow(G),
                        allelesCo = Xco,
                        allelesCa = Xca,
                        OR = NA,
                        Lower.CI = NA,
                        Upper.CI = NA,
                        P = NA,
                        chi.LRT = NA,
                        P.LRT = NA)}
    #} else {res <- NULL}
    #return
    res
  }))


# Output ------------------------------------------------------------------
write.table(output, fileOutput, row.names = FALSE, quote = FALSE)

