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

FROM continuumio/miniconda3:4.7.12

MAINTAINER Tobias Neumann <tobias.neumann.at@gmail.com>

ARG VERSION_ARG

COPY environment.yml /tmp/environment.yml

RUN apt-get update \
    && apt-get install -y procps \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/* \
    && conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda env create --name slamdunk -f /tmp/environment.yml \
    && /opt/conda/envs/slamdunk/bin/pip install git+https://github.com/t-neumann/slamdunk.git@${VERSION_ARG} \
    && rm -rf /opt/conda/pkgs/*

ENV PATH /opt/conda/envs/slamdunk/bin:$PATH
