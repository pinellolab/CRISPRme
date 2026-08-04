# CRISPRme 2.2.0 — Python 3.11 image, built FROM SOURCE.
#
# crispritz 2.8.0 is compiled from source (it is not yet on Bioconda for
# Python 3.11) and crisprme is installed from this build context. The dependency
# pins mirror the from-scratch Python-3.11 validation on ml007 (see PR #131):
#   - azimuth/CRISTA scoring stack: scikit-learn 1.1.3 / numpy 1.24.4 /
#     pandas 2.0.3 / scipy 1.10.1 (the vendored models only unpickle on this combo)
#   - matplotlib-base < 3.9  (matplotlib >= 3.9 runtime-requires numpy >= 1.25,
#     which conflicts with the pinned numpy 1.24.4)
#   - Dash 2.x web stack (dash >= 2.14 bundles the old dash-core/html/renderer/
#     table sub-packages, so those are intentionally dropped)
FROM mambaorg/micromamba

LABEL org.opencontainers.image.authors="ManuelTgn, lucapinello"

ARG crispritz_ref=v2.8.0
ENV SHELL=bash
ENV PREFIX=/opt/conda
USER root

# System build deps to compile the crispritz C++ binaries (g++ + OpenMP + zlib).
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        g++ make git unzip gsl-bin libgsl0-dev libgomp1 zlib1g-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Python 3.11 environment: scoring stack + Dash 2.x web stack + pipeline tools.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        scikit-learn=1.1.3 numpy=1.24.4 scipy=1.10.1 pandas=2.0.3 \
        "matplotlib-base<3.9" \
        biopython more-itertools statsmodels intervaltree \
        pysam bcftools bedtools bedops samtools htslib axel gdown zip gsl \
        importlib-metadata \
        "dash>=2.14,<3" dash-bootstrap-components dash-daq \
        flask flask-caching flask-compress gunicorn werkzeug \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# ---- Build crispritz 2.8.0 from source and install into the env ------------
RUN git clone --depth 1 --branch ${crispritz_ref} \
        https://github.com/pinellolab/CRISPRitz.git /opt/crispritz-src \
    && cd /opt/crispritz-src \
    && g++ -std=c++11 -O3 -fopenmp \
        sourceCode/CRISPR-Cas-Tree/mainParallel.cpp -o buildTST \
    && g++ -std=c++11 -O3 -fopenmp \
        sourceCode/CRISPR-Cas-Tree/searchOnTST.cpp \
        sourceCode/CRISPR-Cas-Tree/detailedOutput.cpp \
        sourceCode/CRISPR-Cas-Tree/convert.cpp \
        -I sourceCode/CRISPR-Cas-Tree/include -o searchTST \
    && g++ -std=c++11 -O3 -fopenmp \
        sourceCode/CRISPRofiler/main.cpp \
        sourceCode/CRISPRofiler/profiling.cpp \
        sourceCode/CRISPRofiler/guide_searching.cpp \
        sourceCode/CRISPRofiler/pam_searching.cpp \
        sourceCode/CRISPRofiler/pre_computation.cpp \
        sourceCode/CRISPRofiler/reading.cpp \
        sourceCode/CRISPRofiler/analysis.cpp -o searchBruteForce \
    && mkdir -p ${PREFIX}/opt/crispritz \
    && cp crispritz.py ${PREFIX}/bin/crispritz.py && chmod +x ${PREFIX}/bin/crispritz.py \
    && cp buildTST searchTST searchBruteForce ${PREFIX}/opt/crispritz/ \
    && cp -R sourceCode/Python_Scripts ${PREFIX}/opt/crispritz/ \
    && rm -rf /opt/crispritz-src

# ---- Install crisprme from this build context (source) ---------------------
# Copy the whole source tree under opt/crisprme/ (mirrors the Bioconda build.sh
# layout) so every module the CLI and web app import is present — PostProcess,
# seq_script, pages, assets, scripts, test, etc. crisprme.py resolves PostProcess
# as <dir-of-crisprme.py>[:-3] + opt/crisprme/, so it must also live in a bin/ dir.
COPY . ${PREFIX}/opt/crisprme/
RUN cp ${PREFIX}/opt/crisprme/crisprme.py ${PREFIX}/bin/crisprme.py \
    && chmod +x ${PREFIX}/bin/crisprme.py \
    && rm -rf ${PREFIX}/opt/crisprme/.git \
    # unzip the CRISTA model at build time (the 276 MB pickle ships zipped in git)
    && if [ -f ${PREFIX}/opt/crisprme/PostProcess/CRISTA_predictors.zip ]; then \
         cd ${PREFIX}/opt/crisprme/PostProcess \
         && unzip -o CRISTA_predictors.zip && rm -f CRISTA_predictors.zip; \
       fi

WORKDIR /root
CMD ["crisprme.py"]
