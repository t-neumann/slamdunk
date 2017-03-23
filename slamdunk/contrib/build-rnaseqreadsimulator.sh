#!/bin/bash

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

wget https://github.com/davidliwei/RNASeqReadSimulator/archive/master.tar.gz -O RNASeqReadSimulator.tar.gz
tar xvzf RNASeqReadSimulator.tar.gz
rm RNASeqReadSimulator.tar.gz

ln RNASeqReadSimulator-master/src/*.py .
