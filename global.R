# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# Workspace ---------------------------------------------------------------
#CRAN packages
# install.packages(c("shiny", "data.table", "dplyr", "lazyeval",
#                    "ggplot2", "ggrepel", "knitr", "markdown","DT",
#                    "UpSetR", "shinythemes", "shinyjs",
#                    "rlang", "scales", "later", "SKAT",
#                    "ggforce", "eulerr", "dndscv", "shiny",
#                    "doParallel"))

library(shiny)
library(shinythemes)
library(shinyjs)
library(data.table)
library(DT)
library(ggplot2)
library(UpSetR)
library(SKAT)
library(ggforce)
library(eulerr)
library(dndscv)
library(doParallel)
library(foreach)

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


# Parallel set up ---------------------------------------------------------
# testing...
#https://www.biostars.org/p/273107/
#https://cran.r-project.org/web/packages/future/
#https://stackoverflow.com/questions/31927035/using-parallel-package-in-shiny
#https://github.com/HenrikBengtsson/doFuture

#cores <- 4
# cores <- makeCluster(detectCores(), type = 'PSOCK')
# 
# if (Sys.info()['sysname'] == 'Windows') {
#   cl <- makeCluster(getOption('cl.cores', cores))
#   registerDoParallel(cl)
#   registerDoSEQ()
#   #on.exit(stopCluster(cl))
# } else {
#   options('mc.cores' = cores)
#   registerDoParallel(cores)
# }


