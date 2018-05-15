####################################################################
# Docker file to build variant2peptide 1.2
# Depends on transvar 2.3.4 & snpeff 4_3t
# This version only support GRCh37/hg19 variant annotation
# Based on ubuntu
####################################################################

# Set the base image to ubuntu
FROM ubuntu:16.04

# File author/Maintainer
MAINTAINER Tenghui Chen

# Set working directory
WORKDIR /opt

# Install the responsible modules
RUN apt-get update \
    && apt-get install -y \
    build-essential \
    ca-certificates \
    python-pip \
    python-dev \
    default-jre \
    default-jdk \
    wget \
    unzip \
    curl \
    git \
    gcc \
    libz-dev \
    zlib1g-dev \
    libbz2-dev \
    libncurses5-dev \
    libncursesw5-dev

## Install the responsible packages
RUN pip install --upgrade pip \
    && pip install \
    configparser \
    future
    
# Install transvar
RUN git clone --recursive https://github.com/zwdzwd/transvar.git \
    && cd transvar \
    && mkdir -p /opt/transvar/lib/python2.7/site-packages \
    && python setup.py install --prefix .

ENV PYTHONPATH $PYTHONPATH:/opt/transvar/lib/python2.7/site-packages/
ENV PATH /opt/transvar/bin:$PATH

# Set up transvar databases
RUN transvar config --download_anno --refversion hg19 \
    && transvar config --download_ref --refversion hg19 \
    && transvar config -k refversion -v hg19

# Setup ENV variables for snpeff
ENV SNPEFF_VERSION=4_3t \
    SNPEFF_HOME=/opt/snpEff

# Install snpEff
RUN wget --quiet -O snpEff_v${SNPEFF_VERSION}_core.zip \
    http://downloads.sourceforge.net/project/snpeff/snpEff_v${SNPEFF_VERSION}_core.zip \
    && unzip snpEff_v${SNPEFF_VERSION}_core.zip -d /opt/ \
    && rm snpEff_v${SNPEFF_VERSION}_core.zip

ENV ANNOT_VERSION=1.2

# Install the annotation script
COPY vcf2peptide_annotation_v${ANNOT_VERSION}.py /usr/local/bin/vcf2peptide_annotation_v${ANNOT_VERSION}.py
RUN chmod 777 /usr/local/bin/vcf2peptide_annotation_v${ANNOT_VERSION}.py

# Change workdir to /usr/local/bin/
WORKDIR /usr/local/bin/
