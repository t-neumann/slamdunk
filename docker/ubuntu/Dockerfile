FROM ubuntu

MAINTAINER Tobias Neumann <tobias.neumann.at@gmail.com>

ENV SLAMDUNK_VERSION 0.1.0
ENV SAMTOOLS_DOWNLOAD_URL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
ENV SLAMDUNK_DOWNLOAD_URL https://github.com/t-neumann/slamdunk.git
ENV PATH /software/slamdunk/bin:$PATH

# binutils is required to run opencl programs
# remove everything that was only needed for building NextGenMap to save image space
RUN buildDeps='git wget gcc g++ libc6-dev make cmake zlib1g-dev python-pip python-dev python-distribute python-pip bzip2 libncurses-dev' \
    runDeps='python default-jre binutils r-base unzip' \
    pipmodules='joblib pysam pybedtools intervaltree pandas numpy biopython' \
    && set -x \
    && apt-get update && apt-get install -y $buildDeps $runDeps --no-install-recommends \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --upgrade pip \
    && pip install $pipmodules \
    && git config --global http.sslVerify false \
    && R -e 'install.packages(c("getopt","ggplot2","gridExtra","RColorBrewer","lattice","matrixStats","dplyr","tidyr"),repos="https://cran.wu.ac.at/")' \
    && mkdir -p /usr/src/ \
    && cd /usr/src/ \
    && wget "$SAMTOOLS_DOWNLOAD_URL" \
    && tar xjvf samtools-1.3.1.tar.bz2 \
    && cd samtools-1.3.1 \
    && make \
    && make prefix=/usr/local/ install \
    && cd .. \
    && git clone "$SLAMDUNK_DOWNLOAD_URL" \
    && cd slamdunk/bin \
    && ./build-ngm.sh \
    && ./build-varscan.sh \
    && ./build-rnaseqreadsimulator.sh \
    && mkdir /software \
    && cp -r /usr/src/slamdunk /software/ \
    && cd /software \
    && rm -r /usr/src/ \
    && apt-get purge -y --auto-remove $buildDeps