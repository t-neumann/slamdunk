#/bin/bash

R_LIBS_SITE=$PWD"/../plot/Rslamdunk"

mkdir -p $R_LIBS_SITE

export R_LIBS_SITE

echo $R_LIBS_SITE

R --no-save <<RSCRIPT

libLoc = .libPaths()[grep("Rslamdunk",.libPaths())]
list.of.packages <- c("getopt","ggplot2","gridExtra","RColorBrewer","lattice","matrixStats","dplyr","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages(lib.loc = libLoc)[,"Package"])]

if(length(new.packages)) install.packages(new.packages, repos="https://cran.wu.ac.at/", lib = libLoc)

RSCRIPT