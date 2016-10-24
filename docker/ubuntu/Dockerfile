FROM ubuntu

MAINTAINER Tobias Neumann <tobias.neumann.at@gmail.com>

ENV SLAMDUNK_VERSION 0.1.0
ENV SLAMDUNK_DOWNLOAD_URL https://github.com/t-neumann/slamdunk.git

# binutils is required to run opencl programs
RUN buildDeps='git wget gcc g++ libc6-dev make cmake zlib1g-dev python-pip python-dev python-distribute python-pip bzip2 libncurses-dev' \
    runDeps='python default-jre binutils r-base unzip' \
    && set -x \
    && apt-get update && apt-get install -y $buildDeps $runDeps --no-install-recommends \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --upgrade pip \
    && cd \
    && git clone $SLAMDUNK_DOWNLOAD_URL \
    && pip install ./slamdunk \
    && rm -R slamdunk \
    && apt-get purge -y --auto-remove $buildDeps