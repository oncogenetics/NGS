# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# Workspace ---------------------------------------------------------------
#CRAN packages
# install.packages(c("shiny", "data.table", "dplyr", "lazyeval",
#                    "ggplot2", "ggrepel", "knitr", "markdown","DT",
#                    "UpSetR", "shinythemes",
#                    "rlang", "scales", "later", "SKAT",
#                    "ggforce", "eulerr"))

library(shiny)
library(data.table)
library(DT)
library(ggplot2)
library(UpSetR)
library(shinythemes)
library(SKAT)
library(ggforce)
library(eulerr)


#Map fonts to Windows
if(Sys.info()['sysname'] == "Windows") {
  windowsFonts(Courier=windowsFont("TT Courier New")) 
}

DToptions <- list(
  pageLength = 50,
  lengthMenu = c(10, 20, 50, 100),
  buttons = list(list(extend = 'csv', filename = 'output')),
  dom = 'Bfrtip',
  # dom = "t",
  fixedHeader = TRUE
  # paging = FALSE
)

# Data --------------------------------------------------------------------
# ANNOT
# GT
# SAMPLE
load("data/data.RData")

# runApp(shinyApp(
#   ui = fluidPage(
#     DT::dataTableOutput("results", height = 300)
#   ),
#   server = function(input, output, session) {
#     output$results <- DT::renderDataTable(
#       mtcars,
#       options = list(pageLength = 200, scrollY = TRUE)
#     )
#   }
# ))
