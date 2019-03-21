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
                   radioButtons("plotParallelVariantGrp", "Parallel plot group:",
                                c("SYMBOL", "IMPACT", "LoF", "REVEL", "CADD_PHRED"),
                                inline = TRUE),
                   splitLayout(
                     sliderInput("qcMissVariant", "Variant missingness", min = 0, max = 0.5, value = 0.5),
                     sliderInput("qcMissSample", "Sample missingness", min = 0, max = 0.5, value = 0.5)
                   ),
                   radioButtons("andOr", "FilterOption", choices = c("AND", "OR"),
                                selected = "OR", inline = TRUE),
                   #MAF
                   sliderInput("ExAC_AF_NFE", "ExAC_AF_NFE *", min = 0, max = 1, value = 0.01),
                   h6("* Non-Finnish European Allele Frequency from ExAC"),
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
                 ), # END input.tabsetPanelMain == 'Variants'
                 conditionalPanel(
                   condition = "input.tabsetPanelMain == 'Sample'",
                   radioButtons("plotParallelSampleGrp", "Parallel plot group:",
                                c("PrCa", "FH", "COD_PrCa", "AgeDiag", "GleasonScore", 
                                  "NCCN", "NICE", "TStage", "NStage", "MStage", "PSADiag"),
                                selected = "AgeDiag",
                                inline = TRUE),
                   checkboxGroupInput("samplePrCa", "Prostate cancer:",
                                      c("0", "1", "NA"),
                                      selected = c("0", "1", "NA"),
                                      inline = TRUE),
                   checkboxGroupInput("sampleEthnicityOA", "Ethnicity OA:",
                                      c("African", "Asian", "European", "Mixed_ethnic", "NA"),
                                      selected = "European",
                                      inline = TRUE),
                   checkboxGroupInput("sampleEthnicity", "Ethnicity Progeny:",
                                      c("0", "1", "6", "7", "8", "NA"),
                                      selected = "1",
                                      inline = TRUE),
                   splitLayout(
                     checkboxGroupInput("sampleFH", "Family history:",
                                        c("0", "1", "NA"), selected = c("0", "1", "NA"), inline = TRUE),
                     checkboxGroupInput("sampleCOD", "COD PrCa:",
                                        c(0, 1, "NA"), selected = c(0, 1, "NA"), inline = TRUE)),
                   sliderInput("sampleAge", "AgeDiag:",
                               min = 0, max = 100, value = c(20, 100), step = 2.5),
                   sliderInput("sampleGleason", "Gleason:",
                               min = 0, max = 10, value = c(2, 10), step = 1),
                   sliderInput("samplePSADiag", "PSA at diagnosis:",
                               min = -1, max = 100, value = c(0, 100), step = 1),
                   checkboxGroupInput("sampleTStage", "TStage:",
                                      c(1:4, "NA"), selected = c(1:4, "NA"), inline = TRUE),
                   splitLayout(
                     checkboxGroupInput("sampleNStage", "NStage:",
                                        c(0, 1, "NA"), selected = c(0, 1, "NA"), inline = TRUE),
                     checkboxGroupInput("sampleMStage", "MStage:",
                                        c(0, 1, "NA"), selected = c(0, 1, "NA"), inline = TRUE))
                 ), # END "input.tabsetPanelMain == 'Sample'"
                 conditionalPanel(
                   condition = "input.tabsetPanelMain == 'Stats'",
                   checkboxGroupInput("cacoType", "Status:",
                                      c("COD_PrCa", "FH", "GleasonScore", "TStage", "NStage", "MStage", "PSADiag", "NCCN", "NICE"),
                                      selected = c("NStage", "MStage"), inline = TRUE)
                 )
                 
                 
                 
                 
    ),#sidebarPanel
    # ~ mainPanel---------------------------------------------------------
    mainPanel(
      tabsetPanel(
        id = "tabsetPanelMain",
        tabPanel("Panel",
                 h4("Panel"),
                 hr(),
                 plotOutput("geneOverlap"), 
                 dataTableOutput("genes")),
        tabPanel("Variants",
                 h4("Variants"),
                 hr(),
                 #tableOutput("testInput"),
                 splitLayout(
                   plotOutput("variantOverlap", width = "80%"),
                   plotOutput("variantSubsetOverlap", width = "80%")),
                 hr(), plotOutput("variantParallel"),
                 hr(), dataTableOutput("annot")),
        tabPanel("Sample",
                 h4("Sample"),
                 hr(),
                 #tableOutput("testFH"),
                 splitLayout(
                   plotOutput("sampleOverlap"),
                   plotOutput("sampleSubsetOverlap")),
                 hr(), plotOutput("phenoParallel"),
                 hr(), plotOutput("phenoNA"),
                 dataTableOutput("pheno")
        ),
        tabPanel("QC",
                 h4("QC"),
                 hr(),
                 tableOutput("testGT"),
                 #plotOutput("qcMissGT"),
                 plotOutput("qcMissGTperGene"),
                 hr(), dataTableOutput("gt")
        ),
        tabPanel("Stats",
                 h4("Stats"),
                 hr(),
                 tabsetPanel(id = "skatBurdenTests", type = "pills",
                             tabPanel("Mutation count", hr(), dataTableOutput("mutCount")),
                             tabPanel("SKAT Gene", hr(), dataTableOutput("geneSkat")),
                             tabPanel("Burden Gene", hr(), dataTableOutput("geneBurden")),
                             tabPanel("SKAT Gene Set", hr(), h3("to-do")),
                             tabPanel("Burden Gene Set", hr(), h3("to-do")))
        ),
        tabPanel("Summary",
                 h4("Summary"),
                 hr(),
                 h3("to-do: some summary plots... ")),
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
