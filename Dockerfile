# Set the base image to micromamba
FROM mambaorg/micromamba

# File Author / Maintainer
LABEL org.opencontainers.image.authors="ManuelTgn, lucapinello"

# Set the variables for version control during installation
ARG crispritz_version=2.7.0
ARG crisprme_version=2.1.12

# set the shell to bash
ENV SHELL bash
# set user as root
USER root

#update packages of the docker image
RUN apt-get update && apt-get install gsl-bin libgsl0-dev -y && apt-get install libgomp1 -y && apt-get clean
RUN apt-get upgrade -y && apt-get clean
RUN apt update
RUN apt upgrade -y

# Install crispritz & crisprme packages
# Python 3.11: the web app now targets the modern Dash stack (Dash 2.x, Flask
# 3.x, Werkzeug 3.x, itsdangerous 2.x). These are intentionally left unpinned
# so conda can resolve a coherent set for Python 3.11.
# NOTE: a full 3.11 image also requires a crispritz Python-3.11 Bioconda build
# (tracked separately); until that lands the crispritz pin here may need to be
# relaxed or sourced from a py3.11 channel.
RUN micromamba install -y -n base -c conda-forge -c bioconda python=3.11 crispritz=$crispritz_version crisprme=$crisprme_version && micromamba clean --all --yes
# Start the base environment
ARG MAMBA_DOCKERFILE_ACTIVATE=1
