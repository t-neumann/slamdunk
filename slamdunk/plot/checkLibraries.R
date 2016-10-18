# Helper function to check whether Rslamdunk libraries are available
# Install if libraries are not available
#
# Author: Tobias Neumann (tobias.neumann.at@gmail.com)
###############################################################################

checkLib <- function(libLoc) {
	
	list.of.packages <- c("getopt","ggplot2","gridExtra","RColorBrewer","lattice","matrixStats","dplyr","tidyr")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages(lib.loc = libLoc)[,"Package"])]
	
	if(length(new.packages)) install.packages(new.packages, repos="http://cran.wu.ac.at/", lib = libLoc)
}