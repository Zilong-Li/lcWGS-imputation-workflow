# Snakemake workflow: lcWGS-imputation-workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)


This workflow is for imputation using low coverage whole genome sequencing data with *QUILT*. Also, it can perform benchmarking for both *QUILT* and *GLIMPSE* given different scenarios.

## Dependencies

- QUILT (QUILT_prepare_reference.R, QUILT.R)
- GLIMPSE v2.0 (GLIMPSE2_split_reference, GLIMPSE2_phase, GLIMPSE2_ligate)
- GLIMPSE v1.1.1 (GLIMPSE_chunk, GLIMPSE_phase, GLIMPSE_ligate)
- samtools
- bcftools

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=Zilong-Li%2FlcWGS-imputation-workflow).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) https://github.com/Zilong-Li/lcWGS-imputation-workflow and its DOI (see above).
