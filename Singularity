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
    # Python 3.11: the web app now targets the modern Dash stack (Dash 2.x,
    # Flask 3.x, Werkzeug 3.x, itsdangerous 2.x). These are intentionally left
    # unpinned so conda can resolve a coherent set for Python 3.11.
    # NOTE: a full 3.11 image also requires a crispritz Python-3.11 Bioconda
    # build (tracked separately); until that lands the crispritz pin here may
    # need to be relaxed or sourced from a py3.11 channel.
    micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        crispritz=2.7.0 \
        crisprme=2.1.12 && \
    micromamba clean --all --yes


%runscript
    echo "Running command in micromamba base..."
    micromamba run -n base "$@"
