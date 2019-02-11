# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# About -------------------------------------------------------------------
# User interface file for shiny


tweaks <- 
  list(
    #push it down 70px to avoid going under navbar
    tags$style(type = "text/css", "body {padding-top: 70px;}"),
    tags$head(
      tags$style(HTML("
                      #graphplotSNP_LDnetwork {
                      border: 1px solid grey;
                      }
                      ")))#,
    #hide red error messages
    # tags$style(type="text/css",
    #            ".shiny-output-error { visibility: hidden; }",
    #            ".shiny-output-error:before { visibility: hidden; }")
      )

# Define UI ---------------------------------------------------------------
shinyUI(
  navbarPage(
    # Application title
    id = "navBarPageID",
    title = div(h4("NGS v0.1",
                   style = "margin-top: 0px;"),
                img(src = "icr_logo_white_on_black.PNG", height = "70px",
                    style = "position: relative; top: -60px; right: -800px;")),
    windowTitle = "NGS",
    fluid = FALSE,
    position = "fixed-top",
    inverse = TRUE,
    
    # 1.Input Data ------------------------------------------------------------
    tabPanel(
      "1.Input Data",
      # ~~~ sidebarPanel ------------------------------------------------------
      sidebarPanel(
        tweaks,
        #Choose data type
        radioButtons("dataType", h4("Input data:"),
                     c("Prostate OncoArray Fine-mapping" = "OncoArrayFineMapping",
                       #"Prostate OncoArray Meta" = "OncoArrayMeta",
                       #"Prostate iCOGS" = "iCOGS",
                       "Custom" = "Custom"
                       #"Example" = "Example"
                     ),
                     selected = "OncoArrayFineMapping")
        ),#sidebarPanel
      # ~~~ mainPanel---------------------------------------------------------
      mainPanel(
        tabsetPanel(
          tabPanel("Summary",
                   h4("Summary"),
                   hr()),
          tabPanel("Association",
                   h4("Association"),
                   hr())
          )#tabsetPanel
        )#mainPanel
      ),#tabPanel - Data
    # 2.Plot Settings ---------------------------------------------------------  
    tabPanel(
      "2.Plot Settings",
      # ~~~ sidebarPanel ----------------------------------------------
             sidebarPanel(
               h4("SNP Filters:")
               ), # END sidebarPanel
      # ~~~ mainPanel -----------------------------------------------
      mainPanel(
        
      ) # END mainPanel
    ), #tabPanel - "2.Plot Settings"
    # 3.Help ------------------------------------------------------------------
    navbarMenu(title = "Help",
               icon = icon("info"),
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
                        includeMarkdown("Markdown/ICR_Data.md"))
               )#navbarMenu - Help
    )#navbarPage
  )#shinyUI



# TESTING -----------------------------------------------------------------
