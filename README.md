# Magma

Bioinformatics pipeline for processing 16S rRNA amplicon sequence data.

## Pipeline

Custom parameters stored in `config.yaml`. See below for a description of the parameters.

## Installation

To use: navigate to project directory. Clone this repository using the following:

```
git clone https://github.com/alanaschick/magma.git magma
```

Note: you need to have **conda** and **snakemake** installed locally in order to run this. See instructions below.

## Executing

Step 1. Enter location of sequences and list of file names in `config.yaml` file.

Step 2. To execute:

Test the pipeline:

```
snakemake -np
```
Run the pipeline:

```
snakemake --use-conda
```

### Remove Adapters

* This pipeline uses the cutadapt software to remove adapters from the raw sequence reads.

### Filtering

* Remove low quality reads using dada2::filterAndTrim function.
* Check quality post filtering with fastqc and multiqc.

### Dada2 Workflow

* Generate error model and denoise reads.
* De-replicate and infer sequence variants.
* Remove bimeras, assign taxonomy.
* Track reads throughout processing, print results to table.

## Conda and Snakemake

To install conda, see the instructions [here](https://github.com/ucvm/synergy/wiki).

To install snakemake using conda, run the following line:

```
conda install -c bioconda -c conda-forge snakemake
```
See the snakemake installation [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for further details.
