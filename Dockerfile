# Set the base image to micromamba
FROM mambaorg/micromamba

# File Author / Maintainer
LABEL org.opencontainers.image.authors="ManuelTgn, lucapinello"

# Set the variables for version control during installation
ARG crispritz_version=2.7.1
ARG crisprme_version=2.1.13

# set the shell to bash
ENV SHELL bash
# set user as root
USER root

#update packages of the docker image
RUN apt-get update && apt-get install gsl-bin libgsl0-dev -y && apt-get install libgomp1 -y && apt-get clean
RUN apt-get upgrade -y && apt-get clean
RUN apt update
RUN apt upgrade -y

# Install crispritz & crisprme packages, plus liftOver (assembly-search's
# liftover dependency -- not yet part of crisprme's own bioconda recipe, so
# listed explicitly here until that recipe is updated)
RUN micromamba install -y -n base -c conda-forge -c bioconda python=3.8 crispritz=$crispritz_version crisprme=$crisprme_version ucsc-liftover && micromamba clean --all --yes
# Start the base environment
ARG MAMBA_DOCKERFILE_ACTIVATE=1
