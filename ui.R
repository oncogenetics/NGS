# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# User interface file for shiny



# Define UI ---------------------------------------------------------------
shinyUI(
  fluidPage(theme = shinytheme("united"),
    titlePanel("NGS v0.1", "NGS"),
    # 1.Input Data ------------------------------------------------------------
    tabPanel(
      "1.Input Data",
      # ~ sidebarPanel ------------------------------------------------------
      sidebarPanel(
        #Choose data type
        uiOutput("ui_data"),
        selectInput("genePanel", "Selected panel(s)", sort(unique(geneListPanel$panel)),
                    selected = "PROCA", multiple = TRUE, selectize = TRUE),
        uiOutput("ui_gene"),
        hr(),
        h5("Filter Variants:"),
        # ClinVar
        selectInput("CLNSIG", "ClinVar", filterCol$CLNSIG,
                    selected = filterCol$CLNSIG[11], multiple = TRUE, selectize = TRUE),
        # VEP
        selectInput("Consequence", "Consequence", filterCol$Consequence,
                    selected = filterCol$Consequence[1], multiple = TRUE, selectize = TRUE),
        checkboxGroupInput("IMPACT", "IMPACT", filterCol$IMPACT, inline = TRUE),
        checkboxGroupInput("LoF", "LoF", filterCol$LoF, inline = TRUE),
        checkboxGroupInput("SIFT", "SIFT", filterCol$SIFT, inline = TRUE),
        checkboxGroupInput("PolyPhen", "PolyPhen", filterCol$PolyPhen, inline = TRUE)
        ),#sidebarPanel
      # ~ mainPanel---------------------------------------------------------
      mainPanel(
        tabsetPanel(
          tabPanel("Variants",
                   h4("Variants"),
                   hr(),
                   dataTableOutput("annot")
          ),
          tabPanel("Association",
                   h4("Association"),
                   hr()),
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
    )#tabPanel
  )#fluidPage
)#shinyUI




# TESTING -----------------------------------------------------------------
