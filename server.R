# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# Server file for shiny

# Define Server -----------------------------------------------------------
shinyServer(function(input, output, session) {

  # geneListVCFselected <- reactive({
  #   sort(unique(geneListVCF[ panel %in% input$data, gene]))
  # })

  subsetGenes <- reactive({
    geneListPanel[ panel %in% input$genePanel, gene]
    # intersect(geneListVCFselected(),
    #           geneListPanel[ panel %in% input$genePanel, gene])
  })



  # ANNOT -------------------------------------------------------------------
  # filter on gene symbol
  subsetAnnot <- reactive({
    # todo: compare MAF over selected VCFs
    # to-do: add gnomadAD link
    # - https://gnomad.broadinstitute.org/variant/1-55516888-G-GA

    rbindlist(
      lapply(input$data, function(i){
        cbind(data = i, ANNOT[[i]][ (SYMBOL %in% input$gene), ])
      }))



    # x <- c("ATM","BRCA1","BRCA2","CHEK2","GEN1","MSH2","MSH5","NBN","PALB2")
    # d <- rbindlist(
    #   lapply(namesVCF[1:2], function(i){
    #     cbind(data = i, ANNOT[[i]][ (SYMBOL %in% x), ])
    #   }))
    # d[ varname %in% d[ duplicated(varname), varname], ][order(varname), ]



  })
  # filter on feature and MAF
  subsetAnnotFilter <- reactive({
    if(length(c(input$CLNSIG, input$Consequence, input$IMPACT,
                input$LoF, input$SIFT, input$PolyPhen)) == 0){
      res <- subsetAnnot()[ ExAC_AF_NFE < input$ExAC_AF_NFE, ]
    } else if(input$andOr == "OR"){
      res <- subsetAnnot()[ ExAC_AF_NFE < input$ExAC_AF_NFE, ][ (
        (grepl(paste(input$CLNSIG, collapse = "|"), CLNSIG) & length(input$CLNSIG) > 0) |
          (grepl(paste(input$Consequence, collapse = "|"), Consequence) & length(input$Consequence) > 0) |
          (IMPACT %in% input$IMPACT & length(input$IMPACT) > 0) |
          (LoF %in% input$LoF & length(input$LoF) > 0) |
          (grepl(paste(input$SIFT, collapse = "|"), SIFT) & length(input$SIFT) > 0) |
          (grepl(paste(input$PolyPhen, collapse = "|"), PolyPhen) & length(input$PolyPhen) > 0)
      ), ]
    } else if(input$andOr == "AND"){
      res <- subsetAnnot()[ ExAC_AF_NFE < input$ExAC_AF_NFE, ][ (
        (grepl(paste(input$CLNSIG, collapse = "|"), CLNSIG) | length(input$CLNSIG) == 0) &
          (grepl(paste(input$Consequence, collapse = "|"), Consequence) | length(input$Consequence) == 0) &
          (IMPACT %in% input$IMPACT | length(input$IMPACT) == 0) &
          (LoF %in% input$LoF | length(input$LoF) == 0) &
          (grepl(paste(input$SIFT, collapse = "|"), SIFT) | length(input$SIFT) == 0) &
          (grepl(paste(input$PolyPhen, collapse = "|"), PolyPhen) | length(input$PolyPhen) == 0)
      ), ]
    }
    # return
    res
  })

  ixVariant <- reactive({
    res <- lapply(ANNOT[ input$data ], function(i){
      #res <- lapply(ANNOT[ namesVCF ], function(i){
      #i="AEPv2"
      merge(i[, list(CHROM, POS, REF, ALT, ix = .I)], subsetAnnotFilter())[, ix]
    })

    res[ lengths(res)!= 0 ]
  })

  # ~Test input --------------------------------------------------------------
  output$testInput <- renderTable({
    t(data.frame(
      nRowSubset = nrow(subsetAnnotFilter()),
      inputCLNSIG = paste(
        paste(input$CLNSIG, collapse = "|"), ":",
        length(input$CLNSIG)),
      inputConsequence = paste(
        paste(input$Consequence, collapse = "|"), ":",
        length(input$Consequence)),
      inputPolyPhen = paste(
        paste(input$PolyPhen, collapse = "|"), ":",
        length(input$PolyPhen))
    ))

  })



  # GT ----------------------------------------------------------------------
  subsetGT <- reactive({
    cols <- c("varname", subsetPheno0()[ , Study.ID])
    dd <- GT[ subsetAnnotFilter()[, varname], ..cols ]

    # missingness on sample and var
    sampleRate <- colSums(is.na(dd[, -1]))/nrow(dd)
    sampleRate <- c("varname", names(sampleRate)[sampleRate < input$qcMissSample])

    variantRate <- which(rowSums(is.na(dd[, -1]))/ncol(dd) < input$qcMissVariant)
    dd[ variantRate, ..sampleRate ]
  })



  # ~Test GT output ----------------------------------------------------------
  output$testGT <- renderTable({
    data.frame(size = paste(dim(subsetGT()), collapse = "x"),
               genes = paste(unique(subsetAnnotFilter()[, SYMBOL]),
                             collapse = ";"))
  })

  # Pheno -------------------------------------------------------------------
  # to-do: get rid of ifelse control... manage NAs better
  subsetPheno0 <- reactive({
    pheno[
      (
        Study.ID %in% unique(unlist(SAMPLE[ input$data ])) &
          ifelse(is.na(PrCa), "NA", as.character(PrCa)) %in% input$samplePrCa &

          (ifelse(is.na(Ethnicity), "NA", Ethnicity)  %in% input$sampleEthnicity |
             ifelse(is.na(Onco_GenoAncestry), "NA", Onco_GenoAncestry)  %in% input$sampleEthnicityOA) &

          ifelse(is.na(FH), "NA", as.character(FH)) %in% input$sampleFH &
          ifelse(is.na(COD_PrCa), "NA", as.character(COD_PrCa)) %in% input$sampleCOD &

          ifelse(is.na(AgeDiag), 0, AgeDiag) >= input$sampleAge[1] &
          ifelse(is.na(AgeDiag), 0, AgeDiag) <= input$sampleAge[2] &

          ifelse(is.na(GleasonScore), 0, GleasonScore) >= input$sampleGleason[1] &
          ifelse(is.na(GleasonScore), 0, GleasonScore) <= input$sampleGleason[2] &

          ifelse(is.na(PSADiag), -1, PSADiag) >= input$samplePSADiag[1] &
          ifelse(ifelse(is.na(PSADiag), -1, PSADiag) > 100, 100,
                 ifelse(is.na(PSADiag), -1, PSADiag)) <= input$samplePSADiag[2] &

          ifelse(is.na(TStage), "NA", as.character(TStage)) %in% input$sampleTStage &
          ifelse(is.na(NStage), "NA", as.character(NStage)) %in% input$sampleNStage &
          ifelse(is.na(MStage), "NA", as.character(MStage)) %in% input$sampleMStage
      ), ]
  })
  subsetPheno <- reactive({
    subsetPheno0()[, list(
      Study.ID,
      Mutation = ifelse(Study.ID %in% studyIDwithMut(), 1, 0),
      PrCa,
      Ethnicity, EthnicityOA,
      AgeDiag, COD_PrCa, FH,
      GleasonScore, Gleason7,
      TStage, NStage, MStage, PSADiag,
      NCCN, NICE) ]
  })

  output$testFH <- renderTable({
    data.frame(class = class(input$sampleFH),
               values = paste(input$sampleFH, collapse = ","))
  })




  # Stats -------------------------------------------------------------------

  # ~Gene SKAT ---------------------------------------------------------------
  #geneSkatConsole

  statSKAT <- reactive({
    d <- merge(subsetAnnotFilter()[ , list(varname = varname, gene = SYMBOL)],
               subsetGT(), by = "varname")

    # d <- merge(subsetAnnot[ , list(varname = varname, gene = SYMBOL)],
    #            subsetGT, by = "varname")
    # setcolorder(d, c("varname", "gene"))

    # out <- foreach(caco = input$cacoType,
    #                      .combine = rbind,
    #                      .export = c("d")
    # ) %dopar% {


    #clusterExport(cl, varlist = "d", envir = as.environment(-1))
    out <- rbindlist(
      #parallel::parLapply(cl, input$cacoType, function(caco){
      lapply(input$cacoType, function(caco){
        # set the trait
        if(caco %in% c("GleasonScore", "TStage", "PSADiag")){
          skatTrait = "C" } else {
            # NStage, MStage
            skatTrait = "D"
          }
        #pp <- subsetPheno[ na.omit(match(tail(colnames(d), -2), Study.ID)), ]
        pp <- subsetPheno()[ na.omit(match(tail(colnames(d), -2), Study.ID)), ]
        ixSample <- which(!is.na(pp[ , ..caco ]))

        #clusterExport(cl, varlist = c("caco"), envir = as.environment(-1))
        rbindlist(
          #parallel::parLapply(cl, split(d, d$gene), function(geneGT){
          lapply(split(d, d$gene), function(geneGT){
            #gene = "BRCA2"; caco = "TStage"
            #geneGT <- split(d, d$gene)[[ gene ]]
            #ixVariant <- which(aa$SYMBOL == gene)

            # SKAT input
            Y <- as.numeric(unlist(pp[ ixSample, ..caco ]))
            if(skatTrait == "D"){ Y <- Y - 1 }
            Z <- t(geneGT[ , -c(1:2) ][ , ..ixSample ])
            Z[ Z == 9 ] <- NA
            X <- pp[ ixSample, AgeDiag]

            obj <- SKAT_Null_Model( Y ~ X, out_type = skatTrait)
            res <- SKAT(Z, obj)

            #output
            if(caco == "PSADiag"){
              x <- table(cut(Y, c(0, 3, 5, Inf)))
            } else {
              x <- table(Y)
            }

            data.table(Gene = head(geneGT$gene, 1),
                       CaCo = caco,
                       CaCoN = paste(paste(names(x), x, sep = "="), collapse = "; "),
                       Trait = skatTrait,
                       P = res$p.value,
                       n.marker = res$param$n.marker,
                       n.marker.test = res$param$n.marker.test)

          }))
      #} #END %dopar%
        })) #END out <- rbindlist(


    #return
    setDT(out)
    out[ order(P), ]


  }) # END statSKAT




  # ~Gene Burden -------------------------------------------------------------
  statBurden <- reactive({
    d <- merge(subsetAnnotFilter()[ , list(varname = varname, gene = SYMBOL)],
               subsetGT(), by = "varname")
    # d <- merge(subsetAnnot[ , list(varname = varname, gene = SYMBOL)],
    #            subsetGT, by = "varname")
    setcolorder(d, c("varname", "gene"))

    d[ is.na(d) ] <- 0

    out <- rbindlist(
      lapply(input$cacoType, function(caco){
        pp <- subsetPheno()[ na.omit(match(tail(colnames(d), -2), Study.ID)), ]
        ixSample <- which(!is.na(pp[ , ..caco ]))

        x <- c(caco, "AgeDiag")
        glmDat <- pp[, ..x]
        colnames(glmDat)[1] <- "caco"
        null.model = glm(caco ~ AgeDiag, family = binomial, data = glmDat)

        rbindlist(
          lapply(split(d, d$gene), function(G){

            gene <- G$gene[1]
            glmDat[, X:=colSums(G[, 3:(ncol(G))] > 0) ]
            Xco = sum(glmDat[ caco == 0, X])
            Xca = sum(glmDat[ caco == 1, X])

            if(sum(glmDat$X, na.rm = TRUE) > 5){
              reg <- glm(caco ~ X + AgeDiag, family = binomial, data = glmDat)
              s.reg <- summary(reg)

              chi.LRT <- 2 * (logLik(reg) - logLik(null.model))
              if(is.na(chi.LRT)){ P.LRT <- NA } else {
                P.LRT <- pchisq(chi.LRT, df = 1, lower.tail = FALSE) }

              Xest <- s.reg$coef["X", "Estimate"]
              Xse <- s.reg$coef["X", "Std. Error"]

              res <- data.frame(gene = gene,
                                caco = caco,
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
              res <- data.frame(gene = gene,
                                caco = caco,
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
      }))

    #return
    out[ order(P), ]
    }) # END statBurden


  # ~dndscv ------------------------------------------------------------
  # http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/dNdScv.html
  dNdScvData <- reactive({
    x <- subsetGT()[ rowSums(subsetGT()[, -1] > 0, na.rm = TRUE) > 5, varname ]
    mutations <- subsetAnnotFilter()[ varname  %in% x, ]

    unique(mutations[, list(
      sampleID = "sample1", #paste0("sample", seq(nrow(mutations))),
      chr = CHROM, pos = POS, ref = REF, mut = ALT) ])
  })

  statdNdScv <- reactive({
    dndsout = dndscv(dNdScvData(),
                     gene_list = input$gene,
                     max_muts_per_gene_per_sample = Inf,
                     max_coding_muts_per_sample = Inf)

    sel_cv <- dndsout$sel_cv
    res <- sel_cv[sel_cv$gene_name %in% input$gene, ]#c("gene_name","qglobal_cv") ]
    rownames(res) <- NULL
    res

    # sel_cv = dndsout$sel_cv
    # print(head(sel_cv), digits = 3)
    #
    # signif_genes = sel_cv[sel_cv$qglobal_cv < 0.1, c("gene_name","qglobal_cv")]
    # rownames(signif_genes) = NULL
    # print(signif_genes)

  })

  # ~Pathway SKAT ------------------------------------------------------------
  statSKATallGenes <- reactive({
    d <- subsetGT()

    out <- rbindlist(
      lapply(input$cacoType, function(caco){
        # set the trait
        if(caco %in% c("GleasonScore", "TStage", "PSADiag")){
          skatTrait = "C" } else {
            # NStage, MStage
            skatTrait = "D"
          }
        pp <- subsetPheno()[ na.omit(match(tail(colnames(d), -2), Study.ID)), ]
        ixSample <- which(!is.na(pp[ , ..caco ]))
        # SKAT input
        Y <- as.numeric(unlist(pp[ ixSample, ..caco ]))
        if(skatTrait == "D"){ Y <- Y - 1 }
        Z <- t(d[ , -1 ][ , ..ixSample ])
        Z[ Z == 9 ] <- NA
        X <- pp[ ixSample, AgeDiag]

        obj <- SKAT_Null_Model( Y ~ X, out_type = skatTrait)
        res <- SKAT(Z, obj)

        #output
        if(caco == "PSADiag"){
          x <- table(cut(Y, c(0, 3, 5, Inf)))
        } else {
          x <- table(Y)
        }

        data.table(Gene = paste0("All genes (n=", length(input$gene), ")"),
                   CaCo = caco,
                   CaCoN = paste(paste(names(x), x, sep = "="), collapse = "; "),
                   Trait = skatTrait,
                   P = res$p.value,
                   n.marker = res$param$n.marker,
                   n.marker.test = res$param$n.marker.test)


      }))

    #return
    out[ order(P), ]


  }) # END statSKATallGenes
  # ~Pathway Burden ----------------------------------------------------------


  # ~Unique mutation count ---------------------------------------------------
  statMutCount <- reactive({
    # subsetGT[, 1:10]
    # subsetAnnot[, 1:10]
    rbindlist(
      lapply(input$cacoType, function(caco){
        #caco = "TStage"
        x <- pheno[ Study.ID %in% colnames(subsetGT()), ]
        res <- colSums(
          sapply(split(x$Study.ID, x[, ..caco ]), function(i){
            (rowSums(subsetGT()[, ..i] > 0, na.rm = TRUE) > 0)*1
          }))
        data.frame(CaCo = caco,
                   CaCoValue = names(res),
                   n = res)
      }))
    })
  studyIDwithMut <- reactive({
    colnames(subsetGT())[ -1 ][ colSums(subsetGT()[, -1]) > 1 ]
  })

  # QC -----------------------------------------------------------------------
  # ~

  # Output DataTables --------------------------------------------------------
  # to-do: rs number link to ensembl:
  # http://www.ensembl.org/Homo_sapiens/Variation/Population?db=core;v=rs4988235
  output$annot <- renderDataTable({
    datatable(
      subsetAnnotFilter(), filter = "top", rownames = FALSE,
      extensions = c("Buttons", "FixedHeader"), options = DToptions)
  })

  output$gt <- renderDataTable({
    d <- subsetGT()[, -1]
    x <- cbind(AA = rowSums(d == 0, na.rm = TRUE),
               AB = rowSums(d == 1, na.rm = TRUE),
               BB = rowSums(d == 2, na.rm = TRUE),
               "NA" = rowSums(is.na(d), na.rm = TRUE))
    datatable(cbind(subsetGT()[, "varname"], x))
  })

  output$pheno <- renderDataTable({
    datatable(
      subsetPheno(), filter = "top", rownames = FALSE,
      extensions = c("Buttons", "FixedHeader"), options = DToptions)
  })

  output$genes <- renderDataTable({
    d <- unique(subsetAnnotFilter()[, .(data, SYMBOL)])
    d <- d[, list(Data = toString(data)), by = SYMBOL ]
    d[ , list(Gene = paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                            SYMBOL,"' target='_blank'>", SYMBOL, "</a>"),
              Data) ]
  }, escape = FALSE)

  output$geneSkat <- renderDataTable({
    datatable(
      statSKAT(),
      filter = "top", rownames = FALSE,
      extensions = c("Buttons", "FixedHeader"), options = DToptions)
  })

  output$allGenesSkat <- renderDataTable({
    datatable(
      statSKATallGenes(),
      filter = "top", rownames = FALSE,
      extensions = c("Buttons", "FixedHeader"), options = DToptions)
  })

  output$geneBurden <- renderDataTable({
    datatable(
      statBurden(),
      #data.frame(x = "testing gene burden output"),
      filter = "top", rownames = FALSE,
      extensions = c("Buttons", "FixedHeader"), options = DToptions)
  })


  output$mutCount <- renderDataTable({
    datatable(
      statMutCount()
      #data.frame(x = "testing output")
      )
  })

  output$dNdScvInput <- renderDataTable({ dNdScvData() })
  output$dNdScv <- renderDataTable({ statdNdScv() })

  #Output Plots -----------------------------------------------------------
  #~GT qc geneSet -------------------------------------------------------
  output$qcMissGT <- renderPlot({
    d <- rbind(
      data.frame(x = "Variants",
                 NArate = rowSums(is.na(subsetGT()[, -1])) / (ncol(subsetGT()) - 1)),
      data.frame(x = "Samples",
                 NArate = colSums(is.na(subsetGT()[, -1])) / (nrow(subsetGT()) - 1)))

    ggplot(d, aes(x = x, y = NArate, col = x)) +
      geom_sina() +
      # geom_violin(trim = FALSE) +
      # geom_jitter(alpha = 0.3) +
      geom_hline(yintercept = input$qcMissVariant, col = "blue", linetype = "dashed" ) +
      geom_hline(yintercept = input$qcMissSample, col = "green", linetype = "dashed" ) +
      scale_y_continuous(name = "Missing Rate", limits = c(0, 1)) +
      ggtitle("Missing rate: Variants and Samples.") +
      theme_classic()
  })
  #~GT qc per gene -----------------------------------------------------
  output$qcMissGTperGene <- renderPlot({
    dd <- merge(subsetAnnot()[ , .(varname, SYMBOL)], subsetGT(), by = "varname")

    d <- rbindlist(lapply(split(dd, dd$SYMBOL), function(i){
      #i=GT[1:1000, 1:20]
      # data.frame(gene = i$SYMBOL[ 1 ],
      #            NArate = sum(is.na(i[, -c("varname", "SYMBOL")]))/prod(dim(i[, -c("varname", "SYMBOL")])))

      cbind(
        gene = i$SYMBOL[ 1 ],
        setNames(data.frame(table(unlist(i[, -c("varname", "SYMBOL")]),
                                  useNA = "always") /
                              prod(dim(i[, -c("varname", "SYMBOL")]))),
                 c("GT", "Freq"))
      )

    }))

    ggplot(d, aes(x = gene, y = Freq, fill = GT)) +
      #geom_boxplot() +
      geom_col() +
      scale_y_continuous(name = NULL, limits = c(0, 1)) +
      ggtitle("Missingness rate per gene") +
      theme_minimal()
  })
  #~Variant Parallel ---------------------------------------------------
  output$variantParallel <- renderPlot({
    # data prep, no Nas, all factor
    d <- subsetAnnotFilter()[ ,list(
      #d <- ANNOT$AEPv2[ ,list(
      varname,
      SYMBOL = factor(SYMBOL),
      CHROM = factor(CHROM, levels = c(1:22, "X", "Y")),
      IMPACT = factor(IMPACT, levels = c("MODIFIER", "LOW", "MODERATE", "HIGH")),
      LoF = factor(LoF, levels = c("LC", "HC")),
      REVEL = cut(REVEL, seq(0, 1, 0.2), labels = c("<0.2", "0.3-0.4", "0.5-0.6", "0.7-0.8", "0.9-1.0")),
      CADD_PHRED = cut(CADD_PHRED, seq(0, 100, 20), labels = c("<20", "30-40", "50-60", "70-80", "90-100")),
      CLNSIG, Consequence
    )]
    # to-do: add: CLNSIG and Consequence
    # d1 <- d[, list(varname, CHROM, SYMBOL, LoF, IMPACT)]
    # Consequence <- d[, lapply(.SD, function(x) unlist(tstrsplit(x, "&", fixed = TRUE))), .SDcols = "Consequence", by = "varname"]
    #
    # CLNSIG <- d[, lapply(.SD, function(x) unlist(tstrsplit(x, ",", fixed = TRUE))), .SDcols = "CLNSIG", by = "varname"
    #             ][, list(varname,  CLNSIG = gsub("^_", "", CLNSIG))]
    #
    # d <- merge(merge(d1, Consequence, by = "varname"), CLNSIG, by = "varname")[ , -1]

    cols <- c("SYMBOL", "IMPACT", "LoF", "REVEL", "CADD_PHRED")
    d[,(cols):= lapply(.SD, function(i){
      factor(ifelse(is.na(i), "NA", as.character(i)),
             levels = c(levels(i), "NA"))
    }), .SDcols = cols]
    d <- d[, ..cols ]

    d <- unique(d[ , value := .N, by = names(d)])

    dd <- gather_set_data(d, seq(ncol(d) - 1))
    dd$x <- factor(dd$x, levels = c("SYMBOL", "IMPACT", "LoF", "REVEL", "CADD_PHRED"))

    ggplot(dd, aes(x, id = id, split = y, value = value)) +
      geom_parallel_sets(aes_string(fill = input$plotParallelVariantGrp),
                         alpha = 0.3, axis.width = 0.3) +
      geom_parallel_sets_axes(axis.width = 0.3) +
      geom_parallel_sets_labels(colour = 'white', angle = 0) +
      scale_x_discrete(name = NULL) +
      ggtitle("Parallel plot: Annotaion relationships",
              paste0("n = ", nrow(subsetAnnotFilter()))
      ) +
      theme_minimal()
  })
  #~Variant overlap ---------------------------------------
  output$variantOverlap <- renderPlot({
    d <- subsetAnnot()[, .(data, varname)]
    d <- split(d$varname, d$data)

    fit <- euler(d, shape = "ellipse")
    plot(fit, quantities = TRUE, main = "Variant overlap")
  })
  output$variantSubsetOverlap <- renderPlot({
    d <- subsetAnnotFilter()[, .(data, varname)]
    d <- split(d$varname, d$data)

    fit <- euler(d, shape = "ellipse")
    plot(fit, quantities = TRUE, main = "Variant subset overlap")
  })


  #~Gene overlap ---------------------------------------
  output$geneOverlap <- renderPlot({
    d <- unique(subsetAnnotFilter()[!is.na(SYMBOL), .(data, SYMBOL)])
    d <- split(d$SYMBOL, d$data)

    fit <- euler(d, shape = "ellipse")
    plot(fit, quantities = TRUE, main = "Gene overlap")
  })

  #~Pheno Parallel -----------------------------------------------------
  output$phenoParallel <- renderPlot({
    # data prep, no Nas, all factor
    #d <- subsetPheno[, list(
    d <- subsetPheno()[, list(
      Mutation = factor(Mutation),
      FH, COD_PrCa,
      AgeDiag = cut(AgeDiag, c(0, 40, 50, 60, 100),
                          labels = c("<40", "40-50", "50-60", ">60")),
      GleasonScore = cut(GleasonScore, c(0, 6, 7, 10),
                         labels = c("1-6", "7", "8-10")),
      Gleason7,
      TStage, NStage, MStage,
      PSADiag = cut(PSADiag, c(0, 3, 5, 20, Inf),
                    labels = c("0-3", "4-5", "6-20", "20+")),
      NCCN, NICE
    )]


    cols <- names(sapply(d, is.factor))[ sapply(d, is.factor) ]
    d[,(cols):= lapply(.SD, function(i){
       factor(ifelse(is.na(as.character(i)), "NA", as.character(i)),
              levels = c(levels(i), "NA"))
     }), .SDcols = cols]

    d <- unique(d[ , value := .N, by = names(d)])

    dd <- gather_set_data(d, seq(ncol(d) - 1))
    dd$x <- factor(dd$x, levels = c(
      "Mutation", "FH", "COD_PrCa", "AgeDiag",
      "GleasonScore", "Gleason7",
      "NCCN", "NICE",
      "TStage", "NStage", "MStage", "PSADiag"))

    ggplot(dd, aes(x, id = id, split = y, value = value)) +
      geom_parallel_sets(aes_string(fill = input$plotParallelSampleGrp),
      #geom_parallel_sets(aes_string(fill = "FH"),
                         alpha = 0.3, axis.width = 0.3) +
      geom_parallel_sets_axes(axis.width = 0.3) +
      geom_parallel_sets_labels(colour = 'white', angle = 0) +
      scale_x_discrete(name = NULL) +
      ggtitle("Parallel plot: Phenotype relationships",
              paste0("n = ", nrow(subsetPheno()))) +
              #paste0("n = ", nrow(subsetPheno))) +
      theme_minimal()



  })
  #~Pheno NAs  -----------------------------------------------------
  output$phenoNA <- renderPlot({
    d <- subsetPheno()[, c("PrCa", "FH", "COD_PrCa", "AgeDiag", "GleasonScore", "NCCN", "NICE", "TStage", "NStage", "MStage", "PSADiag")]
    d <- melt(ifelse(!is.na(d), "Yes", "No"), value.name = "Complete")

    ggplot(d, aes(x = Var2, fill = Complete)) +
      geom_bar() +
      scale_x_discrete(name = NULL) + scale_y_continuous(name = NULL) +
      coord_flip() +
      ggtitle("Penotype data: completeness") +
      theme_minimal()
  })

  #~Sample overlap -------------------------------------------
  output$sampleOverlap <- renderPlot({

    d <- SAMPLE[ input$data ]
    #d <- SAMPLE[ c("AEPv2", "DRG_2441") ]
    fit <- euler(d, shape = "ellipse")
    plot(fit, quantities = TRUE, main = "Sample overlap")
  })

  output$sampleSubsetOverlap <- renderPlot({

    d <- lapply(SAMPLE[ input$data ], function(i){intersect(i, subsetPheno()$Study.ID)})
    #d <- lapply(SAMPLE[ namesVCF[1:2] ], function(i){intersect(i, subsetPheno$Study.ID)})
    fit <- euler(d, shape = "ellipse")
    plot(fit, quantities = TRUE, main = "Sample subset overlap")
  })

  # Dynamic UI --------------------------------------------------------------
  output$ui_genePanel <- renderUI(
    selectInput("genePanel", "Selected panel(s)",
                choices = sort(unique(geneListPanel[ panelType == input$typePanel, panel ])),
                selected = ifelse(input$typePanel == "ICR", "Zsofia_ICR",
                                  sort(unique(geneListPanel[ panelType == input$typePanel, panel ]))[ 1 ]),
                multiple = TRUE, selectize = TRUE))

  output$ui_gene <- renderUI(
    selectInput("gene", "Selected gene(s)", subsetGenes(),
                selected = subsetGenes(), multiple = TRUE, selectize = TRUE))



})#END shinyServer



# TESTING -----------------------------------------------------------------

# testing hilbert ---------------------------------------------------------

#better use it for NGS app, for "fingerprinting mutations"

#library( HilbertVis )

# dataVec <- glmDat$T1-1
# plotLongVector( dataVec, shrinkLength = 100 )
# hMat <- hilbertImage( dataVec,  level = 4)
# showHilbertImage( hMat )
#
# showHilbertImage(rle(dataVec))

# ggplot(points) +
#   geom_voronoi(aes(x,y,fill=distance))
#
#
# ggplot(glmDat, aes(AgeDiag, NCCN, fill = factor(GroupGleasonScore2))) +
#   geom_voronoi_tile() +
#   #geom_voronoi_segment() +
#   geom_point()

