# CRISPRme CI/CD Benchmark Pipeline

This folder provides a continuous integration / continuous deployment (CI/CD) test pipeline to automatically validate the correctness of CRISPRme in retrieving variant-aware off-target sites.

The pipeline is designed to ensure that changes to the codebase do not break correctness, especially with respect to handling population-level genetic variation.

## Purpose of the Pipeline

The CI/CD pipeline verifies that CRISPRme retrieves all off-target sites introduced by genetic variants by comparing its output against a ground-truth benchmark.

### Benchmark strategy

- **Reference genome**: human hg38

- **Variants**: 1000 Genomes Project (1000G)

- **Ground truth**: off-target sites obtained using a **brute-force pairwise sequence alignment [approach](https://github.com/benjaminvyshedskiy/Dynamic_checker)

- **Validation goal**: confirm that CRISPRme retrieves the same set of off-target sites as the brute-force method

The brute-force alignment exhaustively scans the variant-enriched genome and therefore serves as a gold-standard reference, albeit computationally expensive. To keep the CI/CD test tractable, these results are precomputed and reused during the pipeline.

### Scope and Parameters

To ensure scalability while maintaining biological relevance, the benchmark is executed with the following constraints:

- Maximum **4 mismatches**

- Up to **1 DNA bulges**

- Up to **1 RNA bulges**

- PAM: **NGG**

- Guide RNA: **sg1617**

These parameters reflect commonly used thresholds in CRISPR off-target analyses and provide a stringent yet feasible validation setting.

### What the Pipeline Does

When executed, the pipeline automatically:

1. Downloads all required reference data and resources

2. Retrieves precomputed brute-force benchmark results

3. Runs CRISPRme on a **1000G-enriched hg38 reference**

4. Compares CRISPRme results against the benchmark

5. Verifies **100% overlap** between the two sets of off-target sites

No manual intervention is required.

### Requirements

- **Conda / Mamba** (recommended)

- Linux environment

- Internet connection (for downloading reference data)

### How to Run the CI/CD Test Pipeline

**1. Create the CRISPRme environment using Mamba**

Using mamba is strongly recommended for faster and more reliable dependency resolution.

```bash
mamba env create -f environment.yml
```

**2. Activate the environment**
```bash
conda activate crisprme
```

**3. Clone the CRISPRme repository**

The CI/CD test is tied to a fixed version of CRISPRme to ensure reproducibility.

```bash
git clone -b v2.1.9 https://github.com/pinellolab/CRISPRme.git  # currently in branch
cd CRISPRme
```

**4. Run the benchmark script**
```bash
bash test/benchmark/benchmark.sh
```


This single command:

- downloads all required datasets

- runs CRISPRme

- performs the benchmark comparison

- reports the final result

### Interpreting the Results

**Test Passed**

If the pipeline completes successfully, you will see:
```text
Test passed!
```

This indicates that:

- RISPRme retrieved all variant-aware off-target sites

- Results are fully consistent with the brute-force benchmark

### Test Failed

If **any error message** is printed to `stderr` or the script exits prematurely, the test is considered **failed**.

Common causes include:

- missing off-target sites

- mismatches with the benchmark

- corrupted or incomplete input data

- environment or dependency issues

In such cases, the pipeline should be investigated before merging or releasing changes.

### Intended Use

This pipeline is intended for:

- CI/CD validation during development

- Testing before releases

- Reproducibility checks for users

It is not intended as a general-purpose off-target nomination workflow.

### Notes on Reproducibility

- All data sources are versioned and automatically downloaded

- The benchmark relies on precomputed brute-force alignments

- Fixed parameters and guide ensure deterministic behavior

## Suggestions, Feedback, and Custom Benchmarks

We welcome feedback and suggestions to improve the CI/CD benchmark pipeline and CRISPRme more broadly.

### Providing Suggestions and Feedback

If you have suggestions for improving the benchmark, encounter unexpected behavior, or have questions about the CI/CD setup, please:

- Open an issue on the CRISPRme GitHub repository [Issue tracker](https://github.com/pinellolab/CRISPRme/issues)

This helps us track, reproduce, and address issues efficiently.

### Computing Custom Benchmark Datasets

Users interested in validating CRISPRme on custom genomes, variant sets, or parameter configurations (e.g., different organisms, population panels, or mismatch/bulge thresholds) may wish to generate their own benchmark datasets.

Because brute-force alignment approaches are computationally intensive and non-trivial to implement correctly, we strongly recommend contacting the developers before attempting large-scale benchmark generation.

For guidance on:

- generating variant-enriched reference genomes

- selecting appropriate parameters

- running or adapting brute-force alignment methods

- ensuring fair and reproducible comparisons with CRISPRme

please reach out via:

- GitHub [Issues](https://github.com/pinellolab/CRISPRme/issues) (preferred for transparency and reproducibility)

- or by contacting the maintainers directly (see repository contact information)

We are happy to provide methodological guidance, point to existing tools, or advise on best practices to ensure that custom benchmarks are biologically meaningful and technically sound.
