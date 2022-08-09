# Base image
FROM ubuntu:22.04

# Metadata
LABEl base.image="ubuntu:22.04"
LABEL version="1"
LABEL software="Image for ncRNA"
LABEL software.version="20220809"
LABEL description="Image for ncRNA"
LABEL website=""
LABEL documentation=""
LABEL license=""
LABEL tags="Genomics"

# Maintainer
MAINTAINER Mark Ziemann <mark.ziemann@gmail.com>

RUN rm /bin/sh && \
  ln /bin/bash /bin/sh

#numaverage numround numsum
RUN \
  apt-get clean all && \
  apt-get update && \
  apt-get upgrade -y

RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata

RUN \
  apt-get install -y \
    curl \
    nano \
    num-utils \
    wget \
    git \
    perl \
    zip \
    pigz \
    pbzip2 \
    unzip \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libtbb2 \
    default-jdk \
    r-base

########################################
# now downloading a bunch of dependancies
# best to do this in the /sw directory
# also prep where the pipeline will run
# /mnt is for users own data
########################################
RUN mkdir sw

########################################
# Get the ncRNA repo
########################################
RUN git clone git@github.com:casbap/ncRNA.git
ENV DIRPATH /ncRNA
WORKDIR $DIRPATH
RUN chmod -R 777 /ncRNA

########################################
# Get R packages based on own script
########################################

RUN Rscript /ncRNA/docker/rpkgs.R

########################################
# set entrypoint
########################################
#ENTRYPOINT [ "/dee2/code/volunteer_pipeline.sh" ]
