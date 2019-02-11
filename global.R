# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

# Workspace ---------------------------------------------------------------
#CRAN packages
# install.packages(c("shiny", "data.table", "dplyr", "lazyeval",
#                    "ggplot2", "ggrepel", "knitr", "markdown","DT",
#                    "rlang", "scales", "later"))

library(shiny)
library(data.table)

#Map fonts to Windows
if(Sys.info()['sysname'] == "Windows") {
  windowsFonts(Courier=windowsFont("TT Courier New")) 
}


# Data --------------------------------------------------------------------
# ANNOT
# GT
# SAMPLE
load("data/data.RData")
