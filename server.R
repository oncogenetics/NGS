# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# Server file for shiny

# Define Server -----------------------------------------------------------
shinyServer(function(input, output, session) {  
  
  geneListVCFselected <- reactive({
    sort(unique(geneListVCF[ panel %in% input$data, gene]))
  })
  
  subsetGenes <- reactive({
    intersect(geneListVCFselected(),
              geneListPanel[ panel %in% input$genePanel, gene])
  })
  
  
  
  # ANNOT -------------------------------------------------------------------
  # filter on gene symbol
  subsetAnnot <- reactive({
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
  subsetGT <- reactive({
    if(length(ixVariant()) == 0){
      NULL 
    } else {
      gtList <- lapply(names(ixVariant()), function(i){
        d <- data.table(GT[[ i ]][ ixVariant()[[ i ]], , drop = FALSE])
        setnames(d, SAMPLE[[ i ]])
        x <- ANNOT[[ i ]][ixVariant()[[ i ]], ]
        d[ , variant:= paste(x$CHROM, x$POS, x$REF, x$ALT, sep = "_") ]
        d[ , gene:= x$SYMBOL ]
        # d[ , variant:= ifelse(!is.na(x$gnomAD_exomes), x$gnomAD_exomes,
        #                       ifelse(!is.na(x$ID), x$ID,
        #                              paste(x$CHROM, x$POS, x$REF, x$ALT, sep = "_"))) ]
        
      })
      
      if(length(gtList) == 1) {
        res <- gtList[[ 1 ]]
      } else {
        res <- Reduce(function(...) merge(..., by = c("variant", "gene"),
                                          all = TRUE, sort = FALSE), gtList)
        res[ is.na(res) ] <- 9
      }
      # if samples and genotypes overlap between studies, 
      # merge and keep with highst GT: max(-9, 0, 1, 2)
      # to-do: maybe split samples into: repeated and unique, to melt smaller data, then cbind.
      if(any(grepl(".", colnames(res), fixed = TRUE))){
        d <- melt(res, id.vars = c("variant", "gene"))
        d[ , sampleID := tstrsplit(variable, ".", keep = 1, fixed = TRUE) ]
        # missing coded as 9, convert to negative to get max.
        d <- d[ , list(maxValue = max(ifelse(value == 9, -9, value))), by = c("sampleID", "variant", "gene")]
        d[, maxValue:= ifelse(maxValue == -9, 9, maxValue)]
        res <- dcast(d, variant + gene ~ sampleID, value.var = "maxValue")
      }
      
      # return: col1=gene, col2=variantName, cols=samples, rows=vars
      # "variant" column is first
      setcolorder(res, c("gene", "variant"))
      res
      #res[ , c("variant", setdiff(colnames(res), "variant")), with = FALSE]
      #res[ , c("variant", sort(setdiff(colnames(res), "variant")))[1:10], with = FALSE]
      #data.frame(colSums(res[ , c("variant", sort(setdiff(colnames(res), "variant"))), with = FALSE][, -1] == 9))
      #data.frame(table(colnames(res)))
    }
    
  })
  
  # ~Test GT output ----------------------------------------------------------
  output$testGT <- renderTable({
    data.frame(size = paste(dim(subsetGT()), collapse = "x"),
               genes = paste(sort(unique(subsetGT()$gene)), collapse = ";"))
  })
  
  # Pheno -------------------------------------------------------------------
  subsetPheno <- reactive({
    pheno[ Study.ID %in% unlist(SAMPLE[ input$data ]),
           list(
             StudyID = Study.ID,
             #Carrier = as.factor(if_else(Study.ID %in% subsetSamples(), "Yes", "No")),
             rs138213197 = as.factor(rs138213197),
             EthnicityOA = as.factor(Onco_GenoAncestry),
             AgeDiag,
             COD_PrCa = factor(ifelse(is.na(Cause.of.death.is.PrCa), "-",
                                      tolower(Cause.of.death.is.PrCa)), 
                               levels = c("no", "yes", "-")),
             FH = factor(ifelse(is.na(FH), "-", FH),
                         levels = c("No", "Yes", "-")),
             GleasonScore,
             TStage = TStage,
             NStage = NStage,
             MStage = MStage,
             PSADiag = round(PSADiag, 1),
             NCCN = factor(
               ifelse(is.na(NCCN), "-", NCCN),
               levels = c("Low", "Intermediate", "High",
                          "VeryHigh", "Metastatic", "-")),
             NICE = factor(
               ifelse(is.na(NICE), "-", NICE),
               levels = c("Low", "Intermediate", "High", "-"))) ]
    })
  
  # Stats -------------------------------------------------------------------
  
  # ~Gene SKAT ---------------------------------------------------------------
  statSKAT <- reactive({
    
    # aa <- ANNOT$AEPv2
    # gg <- GT$AEPv2
    # ss <- SAMPLE$AEPv2
    #GleasonScore TStage NStage MStage PSADia
    out <- rbindlist(
      lapply(c("GleasonScore", "TStage", "NStage", "MStage", "PSADiag"), function(caco){
      #lapply(c("TStage","NStage", "MStage"), function(caco){
        # set the trait
        if(caco %in% c("GleasonScore", "TStage", "PSADiag")){
          skatTrait = "C" } else {
            # NStage, MStage
            skatTrait = "D"
          }
        pp <- subsetPheno()[ na.omit(match(tail(colnames(subsetGT()), -2), StudyID)), ]
        ixSample <- which(!is.na(pp[ , ..caco ]))
        
        rbindlist(
          lapply(split(subsetGT(), subsetGT()$gene), function(geneGT){
            #gene = "SAMD11"; caco = "GleasonScore"
            #ixVariant <- which(aa$SYMBOL == gene)
            
            # SKAT input
            Y <- unlist(pp[ ixSample, ..caco ])
            Z <- t(geneGT[ , -c(1:2) ][ , ..ixSample ])
            Z[ Z == 9 ] <- NA
            X <- pp[ ixSample, AgeDiag]
            
            obj <- SKAT_Null_Model( Y ~ X, out_type = skatTrait)
            res <- SKAT(Z, obj)
            
            #output
            # data.table(Gene = unique(geneGT$gene),
            #            CaCo = caco,
            #            Trait = skatTrait,
            #            ixSample = length(ixSample),
            #            Y =length(Y),
            #            Z = paste(dim(Z), collapse = "x"),
            #            X = length(X)#paste(dim(X), collapse = "x")
            #            )
            data.table(Gene = unique(geneGT$gene),
                       CaCo = caco,
                       CaCoN = paste(paste(names(table(Y)), table(Y), sep = "="), collapse = "; "),
                       Trait = skatTrait,
                       P = res$p.value,
                       n.marker = res$param$n.marker,
                       n.marker.test = res$param$n.marker.test)
            
          }))
      }))
    
    #return
    out
    
    # testing START ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # aa <- ANNOT$AEPv2
    # gg <- GT$AEPv2
    # ss <- SAMPLE$AEPv2
    # #GleasonScore TStage NStage MStage PSADia
    # out <- rbindlist(
    #   lapply(c("GleasonScore", "TStage", "NStage", "MStage", "PSADiag"), function(caco){
    #     # set the trait
    #     if(caco %in% c("GleasonScore", "TStage", "PSADiag")){
    #       skatTrait = "C" } else {
    #         # NStage, MStage
    #         skatTrait = "D"
    #       }
    #     rbindlist(
    #       lapply(c("SAMD11", "NOC2L", "KLHL17", "PLEKHN1"), function(gene){
    #         #gene = "SAMD11"; caco = "GleasonScore"
    #         ixVariant <- which(aa$SYMBOL == gene)
    #         pp <- pheno[ match(ss, Study.ID), ]
    #         ixSample <- which(!is.na(pp[ , ..caco ]))
    #         
    #         # SKAT input
    #         Y <- unlist(pp[ ixSample, ..caco ])
    #         Z <- t(gg[ ixVariant, ixSample])
    #         Z[ Z == 9 ] <- NA
    #         X <- pp[ ixSample, AgeDiag]
    # 
    #         obj <- SKAT_Null_Model( Y ~ X, out_type = skatTrait)
    #         res <- SKAT(Z, obj)
    #         
    #         #output
    #         data.table(Gene = gene, 
    #                    CaCo = caco,
    #                    Trait = skatTrait,
    #                    P = res$p.value,
    #                    n.marker = res$param$n.marker,
    #                    n.marker.test = res$param$n.marker.test)
    # 
    #       }))
    #   }))
    # testing END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
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
    datatable(subsetGT())
  })
  
  output$pheno <- renderDataTable({
    datatable(
      subsetPheno(), filter = "top", rownames = FALSE,
      extensions = c("Buttons", "FixedHeader"), options = DToptions)
  })
  
  output$skat <- renderDataTable({
    datatable(
      statSKAT(), filter = "top", rownames = FALSE,
      extensions = c("Buttons", "FixedHeader"), options = DToptions)
  })
  # Dynamic UI --------------------------------------------------------------  
  output$ui_genePanel <- renderUI(
    selectInput("genePanel", "Selected panel(s)",
                choices = sort(unique(geneListPanel[ panelType == input$typePanel, panel ])),
                selected = sort(unique(geneListPanel[ panelType == input$typePanel, panel ]))[ 1 ],
                multiple = TRUE, selectize = TRUE))
  
  output$ui_gene <- renderUI(
    selectInput("gene", "Selected gene(s)", geneListVCFselected(),
                selected = subsetGenes(), multiple = TRUE, selectize = TRUE))
  
  
  
})#END shinyServer



# TESTING -----------------------------------------------------------------
