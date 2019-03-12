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
    # todo: ExAC_AF_NFE compare, compare MAF over selected VCFs
    
    rbindlist(
      lapply(input$data, function(i){
        cbind(data = i, ANNOT[[i]][ (SYMBOL %in% input$gene), ])
        
      }))  
  })
  # filter on feature
  subsetAnnotFilter <- reactive({
    if(length(c(input$CLNSIG, input$Consequence, input$IMPACT,
                input$LoF, input$SIFT, input$PolyPhen)) == 0){
      subsetAnnot()  
    } else if(input$andOr == "OR"){
      subsetAnnot()[ (
        (grepl(paste(input$CLNSIG, collapse = "|"), CLNSIG) & length(input$CLNSIG) > 0) |
          (grepl(paste(input$Consequence, collapse = "|"), Consequence) & length(input$Consequence) > 0) |
          (IMPACT %in% input$IMPACT & length(input$IMPACT) > 0) |
          (LoF %in% input$LoF & length(input$LoF) > 0) |
          (grepl(paste(input$SIFT, collapse = "|"), SIFT) & length(input$SIFT) > 0) |
          (grepl(paste(input$PolyPhen, collapse = "|"), PolyPhen) & length(input$PolyPhen) > 0)
      ), ]
    } else if(input$andOr == "AND"){
      subsetAnnot()[ (
        (grepl(paste(input$CLNSIG, collapse = "|"), CLNSIG) | length(input$CLNSIG) == 0) &
          (grepl(paste(input$Consequence, collapse = "|"), Consequence) | length(input$Consequence) == 0) &
          (IMPACT %in% input$IMPACT | length(input$IMPACT) == 0) &
          (LoF %in% input$LoF | length(input$LoF) == 0) &
          (grepl(paste(input$SIFT, collapse = "|"), SIFT) | length(input$SIFT) == 0) &
          (grepl(paste(input$PolyPhen, collapse = "|"), PolyPhen) | length(input$PolyPhen) == 0)
      ), ]
    }
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
  # to-do: sort by chr:pos
  # to-do: remove monomorphic
  # to-do: impute? mean? ref?
  subsetGT <- reactive({
    cols <- c("varname", subsetPheno()[ , StudyID])
    GT[ subsetAnnotFilter()[, varname], ..cols ]
  })
  
  subsetGTmiss <- reactive({
    dd <- subsetGT()
    sampleRate <- colSums(is.na(dd[, -1]))/nrow(dd) 
    sampleRate <- c("varname", names(sampleRate)[sampleRate < input$qcMissSample])
    
    variantRate <- which(rowSums(is.na(dd[, -1]))/ncol(dd) < input$qcMissVariant)
    dd[ variantRate, ..sampleRate ]
  })
  
  
  
  # if(length(ixVariant()) == 0){
  #   NULL 
  # } else {
  #   gtList <- lapply(names(ixVariant()), function(i){
  #     d <- data.table(GT[[ i ]][ ixVariant()[[ i ]], , drop = FALSE])
  #     setnames(d, SAMPLE[[ i ]])
  #     x <- ANNOT[[ i ]][ixVariant()[[ i ]], ]
  #     d[ , variant:= paste(x$CHROM, x$POS, x$REF, x$ALT, sep = "_") ]
  #     d[ , gene:= x$SYMBOL ]
  #     # d[ , variant:= ifelse(!is.na(x$gnomAD_exomes), x$gnomAD_exomes,
  #     #                       ifelse(!is.na(x$ID), x$ID,
  #     #                              paste(x$CHROM, x$POS, x$REF, x$ALT, sep = "_"))) ]
  #     
  #   })
  #   
  #   if(length(gtList) == 1) {
  #     res <- gtList[[ 1 ]]
  #   } else {
  #     res <- Reduce(function(...) merge(..., by = c("variant", "gene"),
  #                                       all = TRUE, sort = FALSE), gtList)
  #     res[ is.na(res) ] <- 9
  #   }
  #   # if samples and genotypes overlap between studies, 
  #   # merge and keep with highst GT: max(-9, 0, 1, 2)
  #   # to-do: maybe split samples into: repeated and unique, to melt smaller data, then cbind.
  #   if(any(grepl(".", colnames(res), fixed = TRUE))){
  #     d <- melt(res, id.vars = c("variant", "gene"))
  #     d[ , sampleID := tstrsplit(variable, ".", keep = 1, fixed = TRUE) ]
  #     # missing coded as 9, convert to negative to get max.
  #     d <- d[ , list(maxValue = max(ifelse(value == 9, -9, value))), by = c("sampleID", "variant", "gene")]
  #     d[, maxValue:= ifelse(maxValue == -9, 9, maxValue)]
  #     res <- dcast(d, variant + gene ~ sampleID, value.var = "maxValue")
  #   }
  #   
  #   # return: col1=gene, col2=variantName, cols=samples, rows=vars
  #   # "variant" column is first
  #   setcolorder(res, c("gene", "variant"))
  #   res
  #   #res[ , c("variant", setdiff(colnames(res), "variant")), with = FALSE]
  #   #res[ , c("variant", sort(setdiff(colnames(res), "variant")))[1:10], with = FALSE]
  #   #data.frame(colSums(res[ , c("variant", sort(setdiff(colnames(res), "variant"))), with = FALSE][, -1] == 9))
  #   #data.frame(table(colnames(res)))
  # }
  
  
  
  # ~Test GT output ----------------------------------------------------------
  output$testGT <- renderTable({
    data.frame(size = paste(dim(subsetGTmiss()), collapse = "x"),
               genes = paste(sort(unique(subsetGTmiss()$gene)), collapse = ";"))
  })
  
  # Pheno -------------------------------------------------------------------
  subsetPheno <- reactive({
    pheno[ Study.ID %in% unique(unlist(SAMPLE[ input$data ])) &
             (Ethnicity %in% input$ethnicity |
                Onco_GenoAncestry %in% input$ethnicityOA),
           list(
             StudyID = Study.ID,
             #Carrier = as.factor(if_else(Study.ID %in% subsetSamples(), "Yes", "No")),
             rs138213197,
             Ethnicity, EthnicityOA,
             AgeDiag, COD_PrCa, FH,
             GleasonScore, TStage, NStage, MStage, PSADiag,
             NCCN, NICE) ]
  })
  
  # Stats -------------------------------------------------------------------
  
  # ~Gene SKAT ---------------------------------------------------------------
  statSKAT <- reactive({
    d <- merge(subsetAnnotFilter()[ , list(varname = varname, gene = SYMBOL)],
               subsetGTmiss(), by = "varname")
    
    # d <- merge(ANNOT$AEPv2[ , list(varname, gene = SYMBOL)],
    #            GT[1:1000, 1:10], by = "varname")
    setcolorder(d, c("varname", "gene"))
    
    out <- rbindlist(
      lapply(c("TStage", "NStage", "MStage"), function(caco){
        #lapply(c("GleasonScore", "TStage", "NStage", "MStage", "PSADiag"), function(caco){
        #caco="TStage"
        #lapply(c("TStage","NStage", "MStage"), function(caco){
        # set the trait
        if(caco %in% c("GleasonScore", "TStage", "PSADiag")){
          skatTrait = "C" } else {
            # NStage, MStage
            skatTrait = "D"
          }
        #pp <- subsetPheno[ na.omit(match(tail(colnames(d), -2), StudyID)), ]
        pp <- subsetPheno()[ na.omit(match(tail(colnames(d), -2), StudyID)), ]
        ixSample <- which(!is.na(pp[ , ..caco ]))
        
        rbindlist(
          lapply(split(d, d$gene), function(geneGT){
            #gene = "LOXL4"; caco = "TStage"
            #geneGT <- split(d, d$gene)[[ gene ]]
            #ixVariant <- which(aa$SYMBOL == gene)
            
            # SKAT input
            Y <- unlist(pp[ ixSample, ..caco ])
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
      }))
    
    #return
    out[ order(P), ]
    
    
  }) # END statSKAT
  
  
  
  # ~Gene Burden -------------------------------------------------------------
  # ~Pathway SKAT ------------------------------------------------------------
  # ~Pathway Burden ----------------------------------------------------------
  
  
  
  # Output ------------------------------------------------------------------
  output$annot <- renderDataTable({
    datatable(
      subsetAnnotFilter(), filter = "top", rownames = FALSE,
      extensions = c("Buttons", "FixedHeader"), options = DToptions)
  })
  
  output$gt <- renderDataTable({
    d <- subsetGTmiss()[, -1]
    x <- cbind(AA = rowSums(d == 0, na.rm = TRUE),
               AB = rowSums(d == 1, na.rm = TRUE),
               BB = rowSums(d == 2, na.rm = TRUE),
               "NA" = rowSums(is.na(d), na.rm = TRUE))
    datatable(cbind(subsetGTmiss()[, "varname"], x))
  })
  
  output$pheno <- renderDataTable({
    datatable(
      subsetPheno(), filter = "top", rownames = FALSE,
      extensions = c("Buttons", "FixedHeader"), options = DToptions)
  })
  
  output$genes <- renderDataTable({
    data.frame(GeneCards = paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                                  subsetGenes(),"' target='_blank'>", subsetGenes(), "</a>"))
  }, escape = FALSE)
  
  # paste0('<a href="', link, gsub("rs", "", id, fixed = TRUE),
  #        '" target="_blank">', id, '</a>')
  
  output$skat <- renderDataTable({
    datatable(
      statSKAT(), filter = "top", rownames = FALSE,
      extensions = c("Buttons", "FixedHeader"), options = DToptions)
  })
  #~Plot: GT qc geneSet -------------------------------------------------------
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
  #~Plot: GT qc per gene -----------------------------------------------------
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
                              prod(dim(i[, -c("varname", "SYMBOL")]))
        ), 
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
  #~Plot: Variant Parallel ---------------------------------------------------
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
      #geom_parallel_sets(aes_string(fill = input$plotParallelGrp), alpha = 0.3, axis.width = 0.1) +
      geom_parallel_sets(aes_string(fill = input$plotParallelVariantGrp), alpha = 0.3, axis.width = 0.1) +
      geom_parallel_sets_axes(axis.width = 0.1) +
      geom_parallel_sets_labels(colour = 'white') +
      scale_x_discrete(name = NULL) +
      ggtitle("Parallel plot: Annotaion relationships", 
              paste0("n = ", nrow(subsetAnnotFilter()))
      ) +
      theme_minimal()
  })
  #~Plot: Variant overlap ---------------------------------------
  output$variantOverlap <- renderPlot({
    d <- subsetAnnotFilter()[, .(data, varname)]
    d <- split(d$varname, d$data)
    
    fit <- euler(d, shape = "ellipse")
    plot(fit, quantities = TRUE)
  })
  
  #~Plot: Gene overlap ---------------------------------------
  output$geneOverlap <- renderPlot({
    d <- unique(subsetAnnotFilter()[!is.na(SYMBOL), .(data, SYMBOL)])
    d <- split(d$SYMBOL, d$data)
    
    fit <- euler(d, shape = "ellipse")
    plot(fit, quantities = TRUE)
  })
  
  #~Plot: Pheno Parallel -----------------------------------------------------
  output$phenoParallel <- renderPlot({
    # data prep, no Nas, all factor
    d <- subsetPheno()[, list(FH, COD_PrCa,
                              AgeDiag = cut(AgeDiag, c(0, 40, 50, 60, 100), 
                                            labels = c("<40", "40-50", "50-60", ">60")),
                              GleasonScore = cut(GleasonScore, c(0, 6, 7, 10),
                                                 labels = c("1-6", "7", "8-10")), 
                              TStage, NStage, MStage,
                              PSADiag = cut(PSADiag, c(0, 3, 5, 20, Inf),
                                            labels = c("0-3", "4-5", "6-20", "20+")),
                              NCCN, NICE
    )]
    
    
    cols <- names(sapply(d, is.factor))[ sapply(d, is.factor) ]
    d[,(cols):= lapply(.SD, function(i){
      factor(ifelse(is.na(i), "NA", as.character(i)),
             levels = c(levels(i), "NA"))
    }), .SDcols = cols]
    
    d <- unique(d[ , value := .N, by = names(d)])
    
    dd <- gather_set_data(d, seq(ncol(d) - 1))
    dd$x <- factor(dd$x, levels = c("FH", "COD_PrCa", "AgeDiag", "GleasonScore", "NCCN", "NICE", "TStage", "NStage", "MStage", "PSADiag"))
    
    ggplot(dd, aes(x, id = id, split = y, value = value)) +
      geom_parallel_sets(aes_string(fill = input$plotParallelSampleGrp), alpha = 0.3, axis.width = 0.1) +
      geom_parallel_sets_axes(axis.width = 0.1) +
      geom_parallel_sets_labels(colour = 'white') +
      scale_x_discrete(name = NULL) +
      ggtitle("Parallel plot: Phenotype relationships", 
              paste0("n = ", nrow(subsetPheno()))) +
      theme_minimal()
  })
  #~Plot: Pheno NAs  -----------------------------------------------------
  output$phenoNA <- renderPlot({
    d <- subsetPheno()[, c("FH", "COD_PrCa", "AgeDiag", "GleasonScore", "NCCN", "NICE", "TStage", "NStage", "MStage", "PSADiag")]
    d <- melt(ifelse(!is.na(d), "Yes", "No"), value.name = "Complete")
    
    ggplot(d, aes(x = Var2, fill = Complete)) +
      geom_bar() +
      scale_x_discrete(name = NULL) + scale_y_continuous(name = NULL) +
      coord_flip() +
      ggtitle("Penotype data: completeness") +
      theme_minimal()
  })
  
  #~Plot: Pheno sample overlap -------------------------------------------
  output$phenoSampleOverlap <- renderPlot({
    
    d <- SAMPLE[ input$data ]
    #d <- SAMPLE[ c("AEPv2", "DRG_2441") ]
    fit <- euler(d, shape = "ellipse")
    plot(fit, quantities = TRUE)
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
