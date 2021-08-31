# archaeodiet

**Pipeline for metagenomic identification of dietary DNA from archaeological samples**.

[![GitHub Actions CI Status](https://github.com/nf-core/archaeodiet/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/archaeodiet/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/archaeodiet/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/archaeodiet/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/archaeodiet.svg)](https://hub.docker.com/r/nfcore/archaeodiet)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible. archaeodiet performs identification and authentication of putative 'dietary' taxa based on DNA from archaeological metagenomic samples.

The pipeline template is based on, and utilises some of the infrastructure of the the [nf-core](https://dx.doi.org/10.1038/s41587-020-0439-x) iniative. Please see the [nf-core website](https://nf-co.re) for more information.


## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility (please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see the [nf-core docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run jfy133/archaeodiet -profile test,<docker/singularity/conda/institute>
```

> This pipeline supports the use of nf-core/configs. Please check the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

iv. Start running your own analysis!

<!-- TODO nf-core: Update the default command above used to run the pipeline -->

```bash
nextflow run jfy133/archaeodiet -profile <docker/singularity/conda/institute> --input '*.fastq.gz'
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The archaeodiet pipeline comes with documentation about the pipeline which you can read at in the [`docs/` directory](docs).

<!-- TODO Add a brief overview of what the pipeline does and how it works -->

## Credits

archaeodiet was originally written by James A. Fellows Yates.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citation

<!-- TODO Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  archaeodiet for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->


## Funding

Development of this pipeline was supported by the ERC Starting Grant project (FoodTransforms) ERC-2015-StG 678901 funded by the European Research Council awarded to Philipp W. Stockhammer (Ludwig Maximilian University, Munich). 
