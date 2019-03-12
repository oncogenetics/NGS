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
                 checkboxGroupInput("data", "Data",
                                    choices = 
                                      setNames(namesVCF, gsub("_", " ", namesVCF)),
                                    selected = namesVCF[1]),
                 hr(),
                 conditionalPanel(
                   condition = "input.tabsetPanelMain == 'Panel'",
                   radioButtons("typePanel", "Select panel type:",
                                c("ICR curated" = "ICR",
                                  "Wood R.D." = "DNA_repair_pathways_Woods",
                                  "MutSigDB Hallmark" = "MutSigDB_hallmark")),
                   uiOutput("ui_genePanel"),
                   uiOutput("ui_gene")
                 ),
                 conditionalPanel(
                   condition = "input.tabsetPanelMain == 'Variants'",
                   splitLayout(
                     sliderInput("qcMissVariant", "Variant missingness", min = 0, max = 0.5, value = 0.5),
                     sliderInput("qcMissSample", "Sample missingness", min = 0, max = 0.5, value = 0.5),
                     radioButtons("andOr", "FilterOption", choices = c("AND", "OR"),
                                  selected = "OR", inline = TRUE)
                   ),
                   
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
                                      inline = TRUE),
                   radioButtons("plotParallelVariantGrp", "Parallel plot group:",
                                c("SYMBOL", "IMPACT", "LoF", "REVEL", "CADD_PHRED"),
                                inline = TRUE)
                 ),
                 conditionalPanel(
                   condition = "input.tabsetPanelMain == 'Phenotype'",
                   radioButtons("plotParallelSampleGrp", "Parallel plot group:",
                                c("FH", "COD_PrCa", "AgeDiag", "GleasonScore", "NCCN", "NICE", "TStage", "NStage", "MStage", "PSADiag"),
                                inline = TRUE),
                   checkboxGroupInput("ethnicityOA", "Ethnicity OA:",
                                      c("African", "Asian", "European", "Mixed_ethnic", NA),
                                      selected = "European",
                                      inline = TRUE),
                   checkboxGroupInput("ethnicity", "Ethnicity:",
                                      c("0", "1", "6", "7", "8", NA),
                                      selected = "1",
                                      inline = TRUE)
                   
                   
                   
                 )
                 
                 
                 
                 
    ),#sidebarPanel
    # ~ mainPanel---------------------------------------------------------
    mainPanel(
      tabsetPanel(
        id = "tabsetPanelMain",
        tabPanel("Panel",
                 h4("Panel"),
                 hr(),
                 dataTableOutput("genes"),
                 plotOutput("geneOverlap")
        ),
        tabPanel("Variants",
                 h4("Variants"),
                 hr(),
                 #tableOutput("testInput"),
                 plotOutput("variantOverlap"),
                 plotOutput("variantParallel"),
                 dataTableOutput("annot")),
        tabPanel("Phenotype",
                 h4("Phenotype"),
                 hr(),
                 plotOutput("phenoSampleOverlap"),
                 plotOutput("phenoParallel"),
                 plotOutput("phenoNA"),
                 dataTableOutput("pheno")
        ),
        tabPanel("Genotype",
                 h4("Genotype"),
                 hr(),
                 tableOutput("testGT"),
                 dataTableOutput("gt"),
                 plotOutput("qcMissGT"),
                 plotOutput("qcMissGTperGene")
                 
        ),
        tabPanel("SKAT",
                 h4("SKAT"),
                 dataTableOutput("skat")),
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
                   tabPanel("Theme",
                            hr("Theme"),
                            shinythemes::themeSelector()),
                   type = "pills")#tabsetPanel
        )#tabPanel
      )#tabsetPanel
    )#mainPanel
  )#fluidPage
)#shinyUI




# TESTING -----------------------------------------------------------------
