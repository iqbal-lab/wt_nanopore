#! /usr/bin/env Rscript
## Script to take the output of alignment_stats and plot coverage across the reference
install.packages('ggplot2')
library(ggplot2)

data <- read.table("~/Dropbox/Nanopore/processed_data/Burn_in/m1/stats/m1_lambda_burn_in_coverage.dat",header=T)
data <- read.table("~/Dropbox/Nanopore/processed_data/Burn_in/m2/stats/m2_lambda_burn_in_coverage.dat",header=T)
data <- read.table("~/Dropbox/Nanopore/processed_data/Burn_in/m3/stats/m3_lambda_burn_in_coverage.dat",header=T)

data <- read.table("~/Dropbox/Nanopore/processed_data/Run_1/Minion1/stats/template_coverage.dat",header=T)
data <- read.table("~/Dropbox/Nanopore/processed_data/Run_1/Minion2/stats/template_coverage.dat",header=T)
data <- read.table("~/Dropbox/Nanopore/processed_data/Run_1/Minion3/stats/template_coverage.dat",header=T)

ggplot(data,aes(x=pos,y=cov)) + geom_line()


data1 <- read.table("~/Dropbox/Nanopore/processed_data/Run_1/Minion1/stats/template_raw_read_stats.dat",header=T)
data2 <- read.table("~/Dropbox/Nanopore/processed_data/Run_1/Minion2/stats/template_raw_read_stats.dat",header=T)
data3 <- read.table("~/Dropbox/Nanopore/processed_data/Run_1/Minion3/stats/template_raw_read_stats.dat",header=T)


ggplot(data1,aes(x=length)) + geom_histogram()
ggplot(data2,aes(x=length)) + geom_histogram()
ggplot(data3,aes(x=length)) + geom_histogram()

