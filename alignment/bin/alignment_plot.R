#! /usr/bin/env Rscript
## Script to take the output of alignment_stats and plot coverage across the reference
# install.packages('ggplot2')
library(ggplot2)
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
dirname <- args[0]
# dirname <- "~/Dropbox/Nanopore/processed_data/Run_1/Minion1/stats/"

printf <- function(...)print(sprintf(...))

coverageData <- read.table(paste(dirname,"template_coverage.dat",sep=""),header=T)
alignmentData <- read.table(paste(dirname,"template_alignment.dat",sep=""),header=T)
## Create output directoy 
outdir <- file.path(dirname, "img")
dir.create(outdir, showWarnings = FALSE)
png(paste(outdir,"position_versus_coverage.png",sep="/"),width=1000,height=800)
ggplot(coverageData,aes(x=pos,y=cov)) + geom_line() + theme_bw() + xlab("Reference Positions") + ylab("Coverage") + ggtitle("Coverage across reference genome")
dev.off()
png(paste(outdir,"histogram_of_read_length.png",sep="/"),width=1000,height=800)
ggplot(alignmentData,aes(x=len)) + geom_histogram(binwidth=100) + theme_bw() + xlab("Read Length") + ggtitle("Histogram of Read Lengths")
dev.off()
## Plot length versus alignment length
png(paste(outdir,"read_length_vs_aligned_read_length.png",sep="/"),width=1000,height=1000)
ggplot(alignmentData,aes(x=len,y=aligned_length)) + geom_point() + theme_bw() + xlab("Read Length") + ylab("Aligned Read Length") + ggtitle("")
dev.off()

numBases <- sum(as.numeric(alignmentData$len))
numAlignedBases <- sum(as.numeric(alignmentData$aligned_length))
M <- sum(as.numeric(alignmentData$M ))
I <- sum(as.numeric(alignmentData$I ))
D <- sum(as.numeric(alignmentData$D ))
S <- sum(as.numeric(alignmentData$S ))
H <- sum(as.numeric(alignmentData$H ))
longestRead <- max(alignmentData$len)
longestAlignedRead <- max(alignmentData$aligned_length)

printf("Number of Bases: %i", numBases)
printf("Number of Aligned Bases : %s", numAlignedBases)
printf("Number of MatchMissMatch Bases : %a", M)
printf("Number of Inserted Bases: %i", I)
printf("Number of Deleted Bases: %i", D)
printf("Number of Softclipped Bases: %i", S)
printf("Number of HardClipped Bases: %i", H)
printf("Longest Read: %i", longestRead)
printf("Longest Aligned Read : %s", longestAlignedRead)
