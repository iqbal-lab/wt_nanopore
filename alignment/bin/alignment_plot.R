#! /usr/bin/env Rscript
## Script to take the output of alignment_stats and plot coverage across the reference
# install.packages('ggplot2')
library(ggplot2)
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
dirname <- args[0]
dirname <- "/Net/wombat/dipro/mmm/data/Nanopore/processed_data/Run_1/Minion1/stats/"
coverageData <- read.table(paste(dirname,"template_coverage.dat",sep=""),header=T)
alignmentData <- read.table(paste(dirname,"template_alignment.dat",sep=""),header=T)
## Create output directoy 
dir.create(file.path(dirname, "img"), showWarnings = FALSE)

ggplot(coverageData,aes(x=pos,y=cov)) + geom_line()
ggplot(alignmentData,aes(x=length)) + geom_histogram()


