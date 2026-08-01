Bootstrap: docker
From: mambaorg/micromamba

%labels
    org.opencontainers.image.authors "ManuelTgn, lucapinello"

%environment
    SHELL=/bin/bash
    # micromamba activate base
    PATH=/opt/conda/bin:$PATH


%post
    echo "Updating system packages..."
    apt-get update && \
    apt-get install -y gsl-bin libgsl0-dev libgomp1 && \
    apt-get upgrade -y && \
    apt-get clean

    echo "Installing CRISPRitz and CRISPRme..."
    micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.8 \
        crispritz=2.7.0 \
        crisprme=2.1.10 && \
    micromamba clean --all --yes


%runscript
    echo "Running command in micromamba base..."
    micromamba run -n base "$@"
