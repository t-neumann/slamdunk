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

FROM ubuntu

MAINTAINER Tobias Neumann <tobias.neumann.at@gmail.com>

ARG VERSION_ARG

# binutils is required to run opencl programs
RUN buildDeps='git wget gcc g++ libc6-dev make cmake zlib1g-dev python-pip python-dev python-distribute python-pip bzip2 libncurses-dev' \
    runDeps='python default-jre binutils r-base unzip' \
    && set -x \
    && apt-get update && apt-get install -y $buildDeps $runDeps --no-install-recommends \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --upgrade pip \
    && pip install git+https://github.com/t-neumann/slamdunk.git@${VERSION_ARG} \
    && apt-get purge -y --auto-remove $buildDeps