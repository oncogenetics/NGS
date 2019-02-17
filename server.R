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
  subsetAnnot <- reactive({
    rbindlist(
      lapply(input$data, function(i){
        cbind(data = i, ANNOT[[i]][ (SYMBOL %in% input$gene), ])
        
      }))  
  })
  
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
  ixVariant <- reactive({
    #     subsetAnnotFilter <- fread("CHROM	POS	REF	ALT
    # 1	203194186	C	T
    #                                1	203194834	C	T
    #                                2	220082505	G	A
    #                                4	89061114	C	T
    #                                5	35037115	C	T
    #                                6	29080004	A	G
    #                                6	29080344	G	A
    #                                7	122635173	A	C
    #                                12	14993439	C	T
    #                                15	28230318	C	T
    #                                X	153713787	C	T
    #                                X	154005148	G	A
    #                                ")
    
    res <- lapply(ANNOT[ input$data ], function(i){
      #res <- lapply(ANNOT[ namesVCF ], function(i){
      #i="AEPv2"
      merge(i[, list(CHROM, POS, REF, ALT, ix = .I)], subsetAnnotFilter())[, ix]
    })
    
    
    res[ lengths(res)!= 0 ]
    #ixVariant <- res
  })
  
  subsetGT <- reactive({
    if(length(ixVariant()) == 0){
      NULL 
    } else {
      #dim(GT[["Eeles_BRCA1_Sanger"]])
      #ANNOT[["Eeles_BRCA1_Sanger"]]
      #dim(GT[["AEPv2"]])
      #subsetAnnotFilter()
      gtList <- lapply(names(ixVariant()), function(i){
        d <- data.table(GT[[ i ]][ ixVariant()[[ i ]], , drop = FALSE])
        setnames(d, SAMPLE[[ i ]])
        x <- ANNOT[[ i ]][ixVariant()[[ i ]], ]
        d[ , variant:= ifelse(!is.na(x$gnomAD_exomes), x$gnomAD_exomes,
                              ifelse(!is.na(x$ID), x$ID,
                                     paste(x$CHROM, x$POS, x$REF, x$ALT, sep = "_"))) ]
        
      })
      
      if(length(gtList) == 1) {
        res <- gtList[[ 1 ]]
      } else {
        res <- Reduce(function(...) merge(..., by = "variant",
                                          all = TRUE, sort = FALSE), gtList)
        res[ is.na(res) ] <- 9
      }
      
      #reverse so that "variant" column is first
      res[ , c("variant", sort(setdiff(colnames(res), "variant")))[1:10], with = FALSE]
    }
    
  })
  
  # ~Test GT output ----------------------------------------------------------
  output$testGT <- renderTable({
    data.frame(dim(subsetGT()))
  })
  
  # Pheno -------------------------------------------------------------------
  subsetPheno <- reactive({
    pheno[ , list(
      StudyID = Study.ID,
      #Carrier = as.factor(if_else(Study.ID %in% subsetSamples(), "Yes", "No")),
      rs138213197 = as.factor(rs138213197),
      EthnicityOA = as.factor(Onco_GenoAncestry),
      AgeDiag,
      COD_PrCa = factor(ifelse(is.na(Cause.of.death.is.PrCa), "-",
                               tolower(Cause.of.death.is.PrCa)), 
                        levels = c("no", "yes", "-")),
      FH = factor(if_else(is.na(FH), "-", FH),
                  levels = c("No", "Yes", "-")),
      GleasonScore,
      TStage = as.factor(TStage),
      NStage = as.factor(NStage),
      MStage = as.factor(MStage),
      PSADiag = round(PSADiag, 1),
      NCCN = factor(
        if_else(is.na(NCCN), "-", NCCN),
        levels = c("Low", "Intermediate", "High",
                   "VeryHigh", "Metastatic", "-")),
      NICE = factor(
        ifelse(is.na(NICE), "-", NICE),
        levels = c("Low", "Intermediate", "High", "-"))) ]
    
  })
  
  
  
  
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
  # Dynamic UI --------------------------------------------------------------  
  output$ui_gene <- renderUI(
    selectInput("gene", "Selected gene(s)", geneListVCFselected(),
                selected = subsetGenes(), multiple = TRUE, selectize = TRUE))
  
  
  
})#END shinyServer



# TESTING -----------------------------------------------------------------
