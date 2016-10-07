#/bin/bash

export R_LIBS_SITE=$PWD"/../plot/Rslamdunk"

echo $R_LIBS_SITE

R --vanilla <<RSCRIPT

libLoc = .libPaths()[grep("Rslamdunk",.libPaths())]
list.of.packages <- c("ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages(lib.loc = libLoc)[,"Package"])]

if(length(new.packages)) install.packages(new.packages, repos="https://cran.wu.ac.at/", lib = libLoc)

RSCRIPT