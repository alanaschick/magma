# Magma

Bioinformatics pipeline for processing 16S rRNA amplicon sequence data.

## Pipeline

Custom parameters stored in `config.yaml`.

### Preprocessing

* Do this thing.
* Do this other thing.

### Filtering

* Remove low quality reads using dada2::filterAndTrim function.
* Check quality post filtering with fastqc and multiqc again.

### Dada2 Workflow

* Generate error model using entire dataset.
* De-replicate and infer sequence variants.
* Remove bimeras, assign taxonomy.
* Track reads throughout processing, print results to table.

