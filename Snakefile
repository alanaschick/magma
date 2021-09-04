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
        "output/multiqc_report.html",
        "output/reads_tracked.txt",
        "output/seqtab_final.rds",
        "output/taxa_final.rds"

rule cutadapt:
    input:
        r1 = config["path"]+"{sample}"+config["for"],
        r2 = config["path"]+"{sample}"+config["rev"]
    params: pre = "{sample}"
    output:
        r1 = "output/temp/cutadapt/{sample}_r1_cutadapt.fastq.gz",
        r2 = "output/temp/cutadapt/{sample}_r2_cutadapt.fastq.gz"
    conda: "utils/envs/cutadapt_env.yaml"
    shell:
            "mkdir -p output/temp/cutadapt/logs; cutadapt -e 0 -O 10 -m 50 -n 2 --discard-untrimmed -g {config[fwd_primer]} "
            "-G {config[rev_primer]} -a {config[rev_primer_rc]} "
            "-A {config[fwd_primer_rc]} -o {output.r1} -p {output.r2} "
            "{input.r1} {input.r2} > output/temp/cutadapt/logs/{params.pre}.cutadapt.log"

rule filter:
    input:
        listfiles = config["list_files"],
        r1 = expand("output/temp/cutadapt/{sample}_r1_cutadapt.fastq.gz", sample=SAMPLES) if config["run_cutadapt"] else expand(config["path"]+"{sample}"+config["for"], sample=SAMPLES),
        r2 = expand("output/temp/cutadapt/{sample}_r2_cutadapt.fastq.gz", sample=SAMPLES) if config["run_cutadapt"] else expand(config["path"]+"{sample}"+config["rev"], sample=SAMPLES)
    output:
        r1 = expand("output/temp/filtered/{sample}_r1_filtered.fastq.gz", sample=SAMPLES),
        r2 = expand("output/temp/filtered/{sample}_r2_filtered.fastq.gz", sample=SAMPLES),
        filt_out = "output/temp/filt_out.rds"
    conda: "utils/envs/dada2_env.yaml"
    script: "utils/scripts/filter.R"

rule fastqc_filt:
    input:
        r1 = "output/temp/filtered/{sample}_r1_filtered.fastq.gz",
        r2 = "output/temp/filtered/{sample}_r2_filtered.fastq.gz"
    output:
        r1 = "output/temp/qc/fastqc/{sample}_r1_filtered_fastqc.html",
        r2 = "output/temp/qc/fastqc/{sample}_r2_filtered_fastqc.html"
    conda: "utils/envs/fastqc_env.yaml"
    shell: "fastqc -o output/temp/qc/fastqc {input.r1} {input.r2}"

rule multiqc:
    input:
        r1 = expand("output/temp/qc/fastqc/{sample}_r1_filtered_fastqc.html", sample=SAMPLES),
        r2 = expand("output/temp/qc/fastqc/{sample}_r2_filtered_fastqc.html", sample=SAMPLES)
    output: "output/multiqc_report.html"
    conda: "utils/envs/multiqc_env.yaml"
    shell: "multiqc -f output/temp -o output -n multiqc_report.html"

rule run_dada2:
    input:
        listfiles = config["list_files"],
        filt_out = "output/temp/filt_out.rds"
    output:
        "output/merged_sequence_lengths.pdf",
        "output/reads_tracked.txt",
        "output/seqtab_final.rds",
        "output/taxa_final.rds"
    conda: "utils/envs/dada2_env.yaml"
    script: "utils/scripts/run_dada2.R"
