#!/usr/bin/env bash
# Install CRISPRme 2.2.0 FROM SOURCE, without the Bioconda crisprme/crispritz
# packages. Builds CRISPRitz 2.8.0 from source and installs both CRISPRitz and
# CRISPRme into the ACTIVE conda environment ($CONDA_PREFIX), using the same
# bin/ + opt/ layout the Bioconda/Docker builds use (crisprme.py resolves
# PostProcess as <dir-of-crisprme.py>[:-3] + opt/crisprme/PostProcess/).
#
# Prerequisites: the crisprme-2.2.0 conda env created from environment.yml and
# ACTIVATED (so a C++ compiler with OpenMP, git and the Python deps are present):
#
#   mamba env create -f environment.yml
#   mamba activate crisprme-2.2.0
#   bash install_from_source.sh
#
# Override the CRISPRitz tag with CRISPRITZ_REF (default v2.8.1).
set -euo pipefail

: "${CONDA_PREFIX:?Activate the conda env first: 'mamba activate crisprme-2.2.0'}"
CRISPRITZ_REF="${CRISPRITZ_REF:-v2.8.1}"
CXX="${CXX:-g++}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # this crisprme checkout

echo ">> Building CRISPRitz ${CRISPRITZ_REF} from source (compiler: ${CXX})"
TMP="$(mktemp -d)"
SRC="${TMP}/crispritz-src"
git clone --depth 1 --branch "${CRISPRITZ_REF}" \
    https://github.com/pinellolab/CRISPRitz.git "${SRC}"
cd "${SRC}"
${CXX} -std=c++11 -O3 -fopenmp \
    sourceCode/CRISPR-Cas-Tree/mainParallel.cpp -o buildTST
${CXX} -std=c++11 -O3 -fopenmp \
    sourceCode/CRISPR-Cas-Tree/searchOnTST.cpp \
    sourceCode/CRISPR-Cas-Tree/detailedOutput.cpp \
    sourceCode/CRISPR-Cas-Tree/convert.cpp \
    -I sourceCode/CRISPR-Cas-Tree/include -o searchTST
${CXX} -std=c++11 -O3 -fopenmp \
    sourceCode/CRISPRofiler/main.cpp \
    sourceCode/CRISPRofiler/profiling.cpp \
    sourceCode/CRISPRofiler/guide_searching.cpp \
    sourceCode/CRISPRofiler/pam_searching.cpp \
    sourceCode/CRISPRofiler/pre_computation.cpp \
    sourceCode/CRISPRofiler/reading.cpp \
    sourceCode/CRISPRofiler/analysis.cpp -o searchBruteForce
# enricher.cpp #includes <zlib.h>; the conda `zlib` pkg installs it under
# $CONDA_PREFIX/include, which is NOT on a system g++'s default path — add it
# explicitly (and the matching lib dir) so the compiled (fast) enricher builds
# regardless of whether $CXX is the conda cxx-compiler or system g++.
${CXX} -std=c++11 -O3 -fopenmp \
    -I"${CONDA_PREFIX}/include" -L"${CONDA_PREFIX}/lib" \
    sourceCode/Python_Scripts/Enrichment/enricher.cpp -o enricher -lz
mkdir -p "${CONDA_PREFIX}/opt/crispritz"
cp crispritz.py "${CONDA_PREFIX}/bin/crispritz.py"
chmod +x "${CONDA_PREFIX}/bin/crispritz.py"
cp buildTST searchTST searchBruteForce "${CONDA_PREFIX}/opt/crispritz/"
cp -R sourceCode/Python_Scripts "${CONDA_PREFIX}/opt/crispritz/"
# Install the compiled binary in both places _enricher_command() (crispritz.py)
# checks, in this order: PATH first, else same-directory as this install
# (corrected_origin_path + "Python_Scripts/Enrichment/enricher"). Without a
# compiled binary in either spot it silently falls back to the pure-Python
# enricher.py -- no error, ~15x slower, see this commit's message for the
# consequence. (add_variants.sh, a separate, not-actually-invoked script, has
# its own unconditional hard-exit if 'enricher' isn't on PATH or beside it --
# same two paths, so satisfied as a side effect regardless.)
cp enricher "${CONDA_PREFIX}/bin/enricher"
chmod +x "${CONDA_PREFIX}/bin/enricher"
cp enricher "${CONDA_PREFIX}/opt/crispritz/Python_Scripts/Enrichment/enricher"
cd "${REPO}"
rm -rf "${TMP}"

echo ">> Installing CRISPRme from ${REPO} into ${CONDA_PREFIX}"
# copy the whole tree so every module the CLI + web app import is present, in
# the layout crisprme.py expects (bin/crisprme.py + opt/crisprme/PostProcess/)
mkdir -p "${CONDA_PREFIX}/opt/crisprme"
cp -R "${REPO}/." "${CONDA_PREFIX}/opt/crisprme/"
rm -rf "${CONDA_PREFIX}/opt/crisprme/.git"
cp "${CONDA_PREFIX}/opt/crisprme/crisprme.py" "${CONDA_PREFIX}/bin/crisprme.py"
chmod +x "${CONDA_PREFIX}/bin/crisprme.py"
# unzip the CRISTA scoring model (ships zipped in git; also auto-unzips on first
# use, but do it now so the first search is not slowed down)
if [ -f "${CONDA_PREFIX}/opt/crisprme/PostProcess/CRISTA_predictors.zip" ]; then
    ( cd "${CONDA_PREFIX}/opt/crisprme/PostProcess" \
        && unzip -o CRISTA_predictors.zip && rm -f CRISTA_predictors.zip )
fi

echo ">> Installed:"
echo "   crisprme.py  -> $(command -v crisprme.py)   ($(crisprme.py --version 2>/dev/null))"
echo "   crispritz.py -> $(command -v crispritz.py)"
echo ">> Done. Try:  crisprme.py complete-test --chrom chr22 --thread 4"
