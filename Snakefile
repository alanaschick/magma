# ******************************
# * Snakefile for 16S pipeline *
# ******************************

# **** Variables ****
configfile: "config.yaml"

# **** Imports ****
import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****
rule all:
    input:
        "results/multiqc_report.html",
        "results/reads_tracked.txt",
        "results/seqtab_final.rds",
        "results/taxa_final.rds"

rule cutadapt:
    input:
        r1 = config["path"]+"{sample}"+config["for"],
        r2 = config["path"]+"{sample}"+config["rev"]
    params: prefix = "{sample}"
    output:
        r1 = "data/trimmed/{sample}_r1_trimmed.fastq.gz",
        r2 = "data/trimmed/{sample}_r2_trimmed.fastq.gz",
        report = "data/qc/cutadapt/{sample}_cutadapt.log"
    conda: "utils/envs/cutadapt_env.yaml"
    shell:
            "cutadapt -e 0 -O 10 -m 50 -n 2 --discard-untrimmed -g {config[fwd_primer]} "
            "-G {config[rev_primer]} -a {config[rev_primer_rc]} "
            "-A {config[fwd_primer_rc]} -o {output.r1} -p {output.r2} "
            "{input.r1} {input.r2} > {output.report}"

rule fastqc_raw:
    input:
        r1 = config["path"]+"{sample}"+config["for"],
        r2 = config["path"]+"{sample}"+config["rev"]
    output:
        r1 = "data/qc/fastqc/{sample}_R1_001_fastqc.html",
        r2 = "data/qc/fastqc/{sample}_R2_001_fastqc.html"
    conda: "utils/envs/fastqc_env.yaml"
    shell: "fastqc -o data/qc/fastqc {input.r1} {input.r2}"

rule multiqc:
    input:
        r1 = expand("data/qc/fastqc/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        r2 = expand("data/qc/fastqc/{sample}_R2_001_fastqc.html", sample=SAMPLES)
    output: "results/multiqc_report.html"
    conda: "utils/envs/multiqc_env.yaml"
    shell: "multiqc -f data/qc -o results -n multiqc_report.html"

rule filter:
    input:
        listfiles = config["list_files"],
        r1 = expand("data/trimmed/{sample}_r1_trimmed.fastq.gz", sample=SAMPLES) if config["run_cutadapt"] else expand(config["path"]+"{sample}"+config["for"], sample=SAMPLES),
        r2 = expand("data/trimmed/{sample}_r2_trimmed.fastq.gz", sample=SAMPLES) if config["run_cutadapt"] else expand(config["path"]+"{sample}"+config["rev"], sample=SAMPLES)
    output:
        r1 = expand("data/filtered/{sample}_r1_filtered.fastq.gz", sample=SAMPLES),
        r2 = expand("data/filtered/{sample}_r2_filtered.fastq.gz", sample=SAMPLES),
        filt_out = "results/filt_out.rds"
    conda: "utils/envs/dada2_env.yaml"
    script: "utils/scripts/filter.R"

rule run_dada2:
    input:
        listfiles = config["list_files"],
        filt_out = "results/filt_out.rds"
    output:
        "results/merged_sequence_lengths.pdf",
        "results/reads_tracked.txt",
        "results/seqtab_final.rds",
        "results/taxa_final.rds"
    conda: "utils/envs/dada2_env.yaml"
    script: "utils/scripts/run_dada2.R"
