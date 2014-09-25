#! /usr/bin/env Rscript
## Script to take the output of raw_read_stats and plot figures
# install.packages('ggplot2')
library(ggplot2)
require(grid)
library(RColorBrewer)

theme2 <- function(){
  return(theme(axis.title= element_text(size = rel(2)),
               axis.text = element_text(size = rel(2)),
               legend.title = element_text(size = rel(2)),
               legend.text = element_text(size = rel(1.5)),
               strip.text = element_text(size = rel(2))) )
}
pal <- c("#E41A1C","#4DAF4A","#377EB8","#984EA3","#FF7F00","#A65628","#F781BF","#999999")
args <- commandArgs(trailingOnly = TRUE)
filename <- args
filename <- "/Users/phelimb/Dropbox/Nanopore/processed_data/lambdaR7/m2/stats/raw_read_stats.dat"
data <- read.table(filename,header=T)
data$GC <- as.numeric(as.character(data$GC))
data$meanQualScore <- as.numeric(as.character(data$meanQualScore))
data$sdQualScore <- as.numeric(as.character(data$sdQualScore))

## Plot a histogram of GC content 
m <- ggplot(data,aes(x=GC)) + geom_histogram(origin = 0) + theme_bw()+ theme2()
m + facet_grid( ~ readType)

ggplot(data,aes(x=length)) + geom_histogram()+ theme_bw()+ theme2() + facet_grid(~ readType)


ggplot(data,aes(x=GC,fill=readType)) + geom_density() + theme_bw()+ theme2()
ggplot(data,aes(x=GC,fill=readType)) + geom_histogram() + theme_bw()+ theme2()

ggplot(data,aes(x=length,y=GC,color=readType)) + geom_point()+ theme_bw()+ theme2()
ggplot(data,aes(x=length,y=meanQualScore)) + geom_point()+ theme_bw()+ theme2()

m <- ggplot(data,aes(x=meanQualScore)) + geom_histogram(binwidth=1) + theme_bw()+ theme2()
m + facet_grid( ~ readType)

ggplot(data,aes(x=meanQualScore,fill=readType)) + geom_density() + theme_bw()+ theme2()

