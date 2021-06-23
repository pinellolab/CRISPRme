# Set the base image to anaconda python 3.8
FROM continuumio/miniconda3

# File Author / Maintainer
# MAINTAINER Samuele Cancelleri

ENV SHELL bash

#update conda channel with bioconda and conda-forge
RUN conda config --add channels defaults 
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

#update packages of the docker system
RUN apt-get update && apt-get install gsl-bin libgsl0-dev -y && apt-get install libgomp1 -y && apt-get clean
RUN apt-get upgrade -y && apt-get clean
RUN apt update
RUN apt upgrade -y

#Install crisprme package
RUN conda install crispritz -y
RUN conda install crisprme -y
RUN conda update crispritz -y
RUN conda update crisprme -y
