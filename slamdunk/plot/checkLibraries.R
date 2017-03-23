# Helper function to check whether Rslamdunk libraries are available
# Install if libraries are not available

# Copyright (c) 2015 Tobias Neumann, Philipp Rescheneder.
#
# This file is part of Slamdunk.
# 
# Slamdunk is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# Slamdunk is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

checkLib <- function(libLoc) {
	
	list.of.packages <- c("getopt","ggplot2","gridExtra","RColorBrewer","lattice","matrixStats","dplyr","tidyr","assertthat","lazyeval","tibble")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages(lib.loc = libLoc)[,"Package"])]
	
	if(length(new.packages)) install.packages(new.packages, repos="http://cran.wu.ac.at/", lib = libLoc, dependencies = TRUE)
}