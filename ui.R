# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# User interface file for shiny



# Define UI ---------------------------------------------------------------
shinyUI(
  fluidPage(
    theme = shinytheme("united"),
    #titlePanel("NGS v0.1", "NGS"),
    titlePanel(title = div("NGS v0.1", img(src = "logoICR.png",
                                           style = "float: right",
                                           height = "40")), "NGS"),
    # ~ sidebarPanel ------------------------------------------------------
    sidebarPanel(width = 3,
      #Choose data type
      checkboxGroupInput("data", "Data",
                         choices = namesVCF, selected = namesVCF[1]),
      conditionalPanel(
        condition = "input.tabsetPanelMain == 'Variants'",
        selectInput("genePanel", "Selected panel(s)",
                    sort(unique(geneListPanel$panel)),
                    selected = "PROCA",
                    multiple = TRUE, selectize = TRUE),
        uiOutput("ui_gene"),
        radioButtons("andOr", "FilterOption", choices = c("AND", "OR"),
                     selected = "OR", inline = TRUE),
        # ClinVar
        selectInput("CLNSIG", "CLNSIG (ClinVar)", filterCol$CLNSIG,
                    #selected = filterCol$CLNSIG[11],
                    multiple = TRUE, selectize = TRUE),
        # VEP
        selectInput("Consequence", "Consequence", filterCol$Consequence,
                    #selected = filterCol$Consequence[1],
                    multiple = TRUE, selectize = TRUE),
        checkboxGroupInput("IMPACT", "IMPACT", filterCol$IMPACT, inline = TRUE),
        checkboxGroupInput("LoF", "LoF", filterCol$LoF, inline = TRUE),
        checkboxGroupInput("SIFT", "SIFT", filterCol$SIFT, inline = TRUE),
        checkboxGroupInput("PolyPhen", "PolyPhen", filterCol$PolyPhen,
                           inline = TRUE)
      )
      
    ),#sidebarPanel
    # ~ mainPanel---------------------------------------------------------
    mainPanel(
      tabsetPanel(
        id = "tabsetPanelMain",
        tabPanel("Variants",
                 h4("Variants"),
                 hr(),
                 tableOutput("testInput"),
                 dataTableOutput("annot")),
        tabPanel("Phenotype",
                 h4("Phenotype"),
                 hr(),
                 dataTableOutput("pheno")),
        tabPanel("Genotype",
                 h4("Genotype"),
                 hr(),
                 tableOutput("testGT"),
                 dataTableOutput("gt")),
        tabPanel("Plots",
                 h4("Plots"),
                 hr()),
        # ~~ Help -------------------------------------------------------
        tabPanel("Help",
                 h4("Help"),
                 hr(),
                 tabsetPanel(
                   tabPanel("About",
                            h4("About"),
                            hr(),
                            includeMarkdown("README.md")),
                   tabPanel("R Session Info",
                            h4("R Session Info"),
                            hr(),
                            includeMarkdown("Markdown/RSessionInfo.md")),
                   tabPanel("Raw data - VCFs",
                            h4("Raw data - VCFs"),
                            hr(),
                            img(src = "logoICR.png"),
                            hr(),
                            includeMarkdown("Markdown/ICR_Data.md")),
                   shinythemes::themeSelector(),
                   type = "pills")#tabsetPanel
        )#tabPanel
      )#tabsetPanel
    )#mainPanel
  )#fluidPage
)#shinyUI




# TESTING -----------------------------------------------------------------
