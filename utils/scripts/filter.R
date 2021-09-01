## Alana Schick
## This is a script to filter 16S microbiome data using dada2's function filterAndTrim
## Note that prior to running this script, primers have been trimmed using Cutadapt - resulting files called samplename_r1_trimmed.fastq.gz
## Last Updated: August 2021

## This is a version of the script to be run using snakemake

library(dada2)
library(tidyverse)

## Set variables
list_of_filenames <- snakemake@input$listfiles

## Set filtering parameters from config file
trimleft <- snakemake@config$trimleft
expected_errors <- c(snakemake@config$expected_errors_forward, snakemake@config$expected_errors_reverse)
truncate <- c(snakemake@config$truncate_forward, snakemake@config$truncate_reverse)

## Cutadapt setting
trimmed <- snakemake@config$run_cutadapt
path_to_raw <- snakemake@config$path

## Make directory for quality plots
dir.create("output/quality_plots/")

#########################################################


## Make a vector of sample names
samples <- scan(list_of_filenames, what = "character")

############# Need to fix this to be dependent on whether or not cutadapt is run
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


# #########################################################
# ######## Step 0: Testing different maxEE values
# results <- array(NA, c(25, 3))
# colnames(results) <- c("for_ee", "rev_ee", "perc")
# count <- 1
# ## Select a sample at random to inspect
# test <- sample(c(1:length(samples)), 1)
#
# for (i in 1:5){
#   for (j in 1:5){
#     out <- filterAndTrim(forward_reads[test],
#                          filtered_forward_reads[test],
#                          reverse_reads[test],
#                          filtered_reverse_reads[test], truncLen=truncate,
#                          maxEE=c(i,j), rm.phix=TRUE,
#                          compress=FALSE, multithread=TRUE, trimLeft = trimleft)
#     res <- out[1,2]/out[1,1]
#     results[count, 1] <- i
#     results[count, 2] <- j
#     results[count,3] <- res
#     count <- count +1
#   }
# }
#
# results <- as.data.frame(results)
#
# gg <- ggplot(results, aes(x = for_ee, y = perc, colour = as.factor(rev_ee))) +
#   geom_point(size = 4) +
#   geom_line(size = 2) +
#   scale_colour_manual(values = rainbow(5, v = 0.8), name = "Error rate Reverse") +
#   xlab("Error rate Forward") +
#   ylab("Percentage reads passed filtering") +
#   ggtitle(paste("Using sample:", samples[test])) +
#   theme_bw()
# gg

#########################################################
####### Step 1: Quality filtering

#out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxEE = expected_errors, multithread = TRUE, rm.phix=TRUE, truncLen=truncate, trimLeft = trimleft, compress = TRUE)

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





