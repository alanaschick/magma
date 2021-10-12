## Alana Schick
## This is a script to process 16S microbiome data using dada2
## Note that prior to running this script, primers have been trimmed using Cutadapt and reads have been filtered
## Last updated: OCTOBER 12 2021

## This is a version of the script to be run using snakemake

library(dada2)
#packageVersion("dada2")
library(tidyverse)

## Set variables
list_of_filenames <- snakemake@input$listfiles
tax_ref <- snakemake@config$dada2_tax_ref
spe_ref <- snakemake@config$dada2_spe_ref
bacterial <- snakemake@config$bacterial

## Make a vector of sample names
samples <- scan(list_of_filenames, what = "character")

out <- readRDS(snakemake@input$filt_out)


## file names of forward and reverse reads, after quality filtering
filtered_forward_reads <- file.path("output/temp/filtered", paste0(samples, "_r1_filtered.fastq.gz"))
filtered_reverse_reads <- file.path("output/temp/filtered", paste0(samples, "_r2_filtered.fastq.gz"))


####### Step 2: Generate error model of data

err_forward <- learnErrors(filtered_forward_reads, multithread = TRUE)
err_reverse <- learnErrors(filtered_reverse_reads, multithread = TRUE)



####### Step 3: Derepliate sequences

derep_forward <- derepFastq(filtered_forward_reads, verbose = TRUE)
derep_reverse <- derepFastq(filtered_reverse_reads, verbose = TRUE)

names(derep_forward) <- samples
names(derep_reverse) <- samples



####### Step 4: Infer sequence variants

dadaF <- dada(derep_forward, err = err_forward, multithread = TRUE)
dadaR <- dada(derep_reverse, err = err_reverse, multithread = TRUE)



####### Step 5: Merge forward and reverse reads
merged <- mergePairs(dadaF, derep_forward, dadaR, derep_reverse, verbose = TRUE)



####### Step 6: Generate count table
seqtab <- makeSequenceTable(merged)




####### Step 7: Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)


## How are the reads concentrated in the merged sequence lengths
readsbyseqlen <- tapply(colSums(seqtab.nochim), nchar(colnames(seqtab.nochim)),sum)
pdf(file.path("output", "merged_sequence_lengths.pdf"))
plot(as.integer(names(readsbyseqlen)),readsbyseqlen, xlab = "Merged length", ylab = "Total reads")
dev.off()


####### Step 7b: Track reads throughout processing

getN <- function(x) sum (getUniques(x))
summary_tab <- data.frame(row.names=samples, Input=out[,1], Filtered=out[,2], Denoised=sapply(dadaF, getN), Merged=sapply(merged, getN), Non.Chimeric=rowSums(seqtab.nochim), Total.Perc.Remaining = round(rowSums(seqtab.nochim)/out[,1]*100,1))


pdf(file.path("output", "reads_remaining.pdf"))
hist(summary_tab$Total.Perc.Remaining, col = "blue", breaks = 50, main = "Total", xlab = "Percentage reads remaining")
dev.off()

## Write this table to output
write.table(summary_tab, file.path("output", "reads_tracked.txt"))

summary_tab$Sample <- rownames(summary_tab) 
summary_tab <- summary_tab %>% separate(Sample, c("Sample", "temp"), sep = "_S") 
summary_tab$Sample <- factor(summary_tab$Sample, levels = summary_tab$Sample[order(summary_tab$Non.Chimeric)])
summary_tab_long <- summary_tab %>% gather("QC.Step", "Reads", Input:Non.Chimeric)
summary_tab_long$QC.Step <- factor(summary_tab_long$QC.Step, levels = c("Input", "Filtered", "Denoised", "Merged", "Non.Chimeric"))


gg <- ggplot(summary_tab_long, aes(x = Sample, y = Reads, color = QC.Step)) +
	geom_point(size = 2) +
	scale_color_manual(values = rainbow(5, v = 0.8)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90))
pdf(file.path("output", "read_tracking.pdf"))
gg
dev.off()



####### Step 8: Assign Taxonomy
################################################
## To do: add parameter to config file to allow user to decide whether or not to allow multiples
################################################
taxa <- assignTaxonomy(seqtab.nochim, tax_ref, multithread = TRUE)

## Add Species (if 16S, not if ITS)
if (bacterial == T){
	taxa <- addSpecies (taxa, spe_ref, allowMultiple = TRUE)
}	


####### Step 9: Save output

saveRDS(seqtab.nochim, file.path("output", "seqtab_final.rds"))
saveRDS(taxa, file.path("output", "taxa_final.rds"))





