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
  
  ixVariant <- reactive({
    lapply(input$data, function(i){
      which(ANNOT[[ i ]][ SYMBOL %in% subsetGenes(), ])
    })
  })
  
  # ANNOT -------------------------------------------------------------------
  subsetAnnot <- reactive({
    rbindlist(
      lapply(input$data, function(i){
        cbind(data = i, ANNOT[[i]][ (SYMBOL %in% subsetGenes()), ])
      
    }))  
  })

  
  
  # GT ----------------------------------------------------------------------
  
  
  # Output ------------------------------------------------------------------

  output$annot <- renderDataTable({ datatable(
    subsetAnnot(), filter = "top", rownames = FALSE,
    extensions = "Buttons", options = DToptions) })
  
  
  # Dynamic UI --------------------------------------------------------------  
  output$ui_data <- renderUI({
    checkboxGroupInput("data", "Data",
                       choices = namesVCF, selected = namesVCF[1])
  })

  output$ui_gene <- renderUI(
    selectInput("gene", "Selected gene(s)", geneListVCFselected(),
                selected = subsetGenes(), multiple = TRUE, selectize = TRUE))
  
  
  
})#END shinyServer



# TESTING -----------------------------------------------------------------
