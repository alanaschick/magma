## Alana Schick
## This is a script to filter 16S microbiome data using dada2's function filterAndTrim
## Note that prior to running this script, primers have been trimmed using Cutadapt - resulting files called samplename_r1_trimmed.fastq.gz
## Last Updated: Nov 22 2021

## This is a version of the script to be run using snakemake

library(dada2)
library(tidyverse)
library(colorRamps)

## Set variables
list_of_filenames <- snakemake@input$listfiles

## Set filtering parameters from config file
trimleft <- c(snakemake@config$trimleft_forward, snakemake@config$trimleft_reverse)
expected_errors <- c(snakemake@config$expected_errors_forward, snakemake@config$expected_errors_reverse)
truncate <- c(snakemake@config$truncate_forward, snakemake@config$truncate_reverse)
readlength <- snakemake@config$readlength

## Cutadapt setting
trimmed <- snakemake@config$run_cutadapt
path_to_raw <- snakemake@config$path

## Exploring parameter space options
trunc_param <- snakemake@config$explore_truncation_params
ee_param <- snakemake@config$explore_quality_params

## Make directory for quality plots
dir.create("output/quality_plots/")

#########################################################


## Make a vector of sample names
samples <- scan(list_of_filenames, what = "character")

## file names of forward and reverse reads, before quality filtering
if (trimmed == T){
	forward_reads <- file.path("output/temp/cutadapt", paste0(samples, "_r1_cutadapt.fastq.gz"))
	reverse_reads <- file.path("output/temp/cutadapt", paste0(samples, "_r2_cutadapt.fastq.gz"))
}

if (trimmed == F){
	forward_reads <- file.path(path_to_raw, paste0(samples, "_R1_001.fastq.gz"))
	reverse_reads <- file.path(path_to_raw, paste0(samples, "_R2_001.fastq.gz"))
}

## file names of forward and reverse reads, after quality filtering
filtered_forward_reads <- file.path("output/temp/filtered", paste0(samples, "_r1_filtered.fastq.gz"))
filtered_reverse_reads <- file.path("output/temp/filtered", paste0(samples, "_r2_filtered.fastq.gz"))

#########################################################
######## Step 0: Exploring filtering parameters

if (trunc_param == T){
## Explore truncation parameters
results <- NULL
## Select a set of samples at random to inspect
test <- sample(c(1:length(samples)), 12)
 
for (i in seq(from = readlength-45, to = readlength, by = 5)){
  for (j in seq(from = readlength-90, to = readlength, by = 10)){
  	truncparam <- c()
    out <- filterAndTrim(forward_reads[test],
                         filtered_forward_reads[test],
                         reverse_reads[test],
                         filtered_reverse_reads[test], 
						 truncLen=c(i,j),
                         maxEE=expected_errors, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, trimLeft = trimleft)
    res <- data.frame(Sample = rownames(out), perc = out[,2]/out[,1], for_trunc = i, rev_trunc = j)
    results <- rbind(results, res)
  }
}

results <- results %>% separate(Sample, c("Name", "Sample"), sep = "_S")

gg <- ggplot(results, aes(x = for_trunc, y = perc, colour = as.factor(rev_trunc))) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_colour_manual(values = sample(primary.colors(20), 10), name = "Truncate Reverse") +
  xlab("Truncate Forward") +
  ylab("Percentage reads passed filtering") +
  ggtitle("Truncation parameters") +
  theme_minimal() +
  facet_wrap(~Name)


pdf(file.path("output", "truncation_parameters.pdf"))
print(gg)
dev.off()
}

if (ee_param == T){
## Explore expected error values
results <- NULL

for (i in 1:5){
	for (j in 1:5){
		out <- filterAndTrim(forward_reads[test], 
							filtered_forward_reads[test], 
							reverse_reads[test],
							filtered_reverse_reads[test], 
							truncLen=truncate,
							maxEE=c(i,j), rm.phix=TRUE,
							compress=TRUE, multithread=TRUE, trimLeft = trimleft)
		res <- data.frame(Sample = rownames(out), perc = out[,2]/out[,1], for_error = i, rev_error = j)
 		results <- rbind(results, res)

	}
}

results <- results %>% separate(Sample, c("Name", "Sample"), sep = "_S")

gg <- ggplot(results, aes(x = for_error, y = perc, colour = as.factor(rev_error))) +
	geom_point(size = 2) +
	geom_line(size = 1) +
	scale_colour_manual(values = rainbow(5, v = 0.8), name = "Error Rate Reverse") +
	xlab("Error Rate Forward") +
	ylab("Percentage reads passed filtering") +
	ggtitle("Expected error parameters") +
	theme_minimal() +
	facet_wrap(~Name)


pdf(file.path("output", "expected_error_parameters.pdf"))
print(gg)
dev.off()
}

#########################################################
####### Step 1: Quality filtering

out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxEE = expected_errors, multithread = TRUE, rm.phix=TRUE, trimLeft = trimleft, compress = TRUE, truncLen = truncate)

saveRDS(out, file.path("output/temp", "filt_out.rds"))

####### Step 2: Plot quality profiles
## Select 10 samples at random to inspect/print
toplot <- sample(c(1:length(samples)), 10)

## Plotting forward versus reverse quality
for (i in 1:10){
  pdf(paste("output/quality_plots/quality_", samples[toplot[i]], ".pdf", sep = ""))
  print(plotQualityProfile(c(filtered_forward_reads[toplot[i]], forward_reads[toplot[i]], filtered_reverse_reads[toplot[i]], reverse_reads[toplot[i]])))
  dev.off()
}





