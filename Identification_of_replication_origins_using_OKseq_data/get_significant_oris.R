#!/usr/bin/Rscript

# Script takes the outputfile .*deltarfs.txt from detect_strandswitches.pl as input.
# Script takes inputfilename from the commandline. 
# Script outputs estimated oris (that had a benjamini hochberg adjusted pvalue < 0.01 )
# Input: *deltarfs.txt from detect_strandswitches.pl
# Output: *estimated_oris.txt contains information on the significant origins
# Outputformat: chromosome   start(=origin start)   stop(=origin stop)  ratio(=delta_fr ratio)  deltafr leftf(=number of fragments from forward strand to left of origin)   rightr(=number of fragments from reverse strand to right of origin)	p-value(from binomialtest)	rank	bha-adjusted p-value. Columns are tab-separated.
#
# Author: Thamar Jessurun Lobo
# Date: 16-01-2019

# Imports
require(data.table)
library(ggplot2)

# Args from commandline
args = commandArgs(trailingOnly=TRUE)
print(args)
print(length(args))

# Function outputs origin locations for significant strandswitches with a positive deltaRF value (a significant negative deltafr value indicates a fork merger zone)
get_significant_oris <- function(x) {
    file = args[x]
    print(file)
    scale = as.integer(sub(".*bam_(\\d{1,})[nr]deltarfs.*", '\\1', basename(file))) # Extract scale from filename
    n_unique_reads =  # number of valid reads as reported by detect_strandswitches.pl
    pcutoff = 0.01 # BH adjusted p-values should be lower then pcutoff. Default is 0.01

    # Read in data
    dt <- fread(file, sep="\t", col.names=c("chrom", "start", "stop", "ratio", "deltafr", "leftf", "rightr"))

    # Binomal test gives p value
    dt$p_value <- sapply(1:nrow(dt), FUN=function(x) {return(binom.test(c(dt[x, (deltafr+scale)], (2*scale) - dt[x, (deltafr+scale)]))$p.value)} )

    # Get rank pvalues
    dt <-dt[order(dt$p_value),]
    dt$rank <- rank(dt$p_value, ties.method = "min")

    # Get Benjamini hochberg adjusted p value
    dt$bha <- dt$p_value * (n_unique_reads/dt$rank)
    head(dt)

    # Select origin switches 
    oris <- dt[which(dt$bha < pcutoff & dt$deltarf > 0),]
    print(paste(sep=" ", "Number of significant oris detected:", nrow(oris)))

    # Sort origins on location
    oris <- oris[order(oris$chrom, oris$start, oris$stop),]

	# Write origins to outputfile
    write.table(oris, file=paste(sep="_", sub(".txt", "", file), "estimated_oris.txt"), sep="\t", quote = F, col.names = F, row.names = F)
}

# For each command line argument, run get_significant_oris
junk <- lapply(1:length(args), FUN=get_significant_oris)
