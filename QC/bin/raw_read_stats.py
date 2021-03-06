#! /usr/bin/env python
import sys
import os
sys.path.append('../src')
sys.path.append('../../fast5_file_manipulation/src')
from fast5 import fast5File
from optparse import OptionParser
import glob
import logging
import csv
from Bio import SeqIO
# from Bio.SeqUtils import GC
import numpy as np
import operator
def GC(seq):
	"""Calcultes GC/ATCG content in a sequence"""
	A = seq.count('A')
	T = seq.count('T')
	C = seq.count('C')
	G = seq.count('G')
	return float(C+G) / float(A+T+G+C)

## Takes a directory of basecalled ONT fast5 files or fastafiles and returns some basic stats 
## on the number of bases, reads, quality scores etc. 
parser = OptionParser()
parser.description = """
Calculates some basic stats on ONT fastq reads"""
parser.add_option("-d", "--dir", dest="directory",
                  help="path to directory containing fast5 or fastq files")
parser.add_option("-o", "--out", dest="outfile_name",
                  help="path to output",default="./")
(opt, args) = parser.parse_args()
if opt.directory is None:
    parser.print_help()
    sys.exit()

### Get all the files in the directory 
fqFilelist = glob.glob(opt.directory+"/*.fq") + glob.glob(opt.directory+"/*.fastq")
fast5Filelist = glob.glob(opt.directory+"/*.fast5") 
if fast5Filelist:
	logging.error("Haven't written fast5 handler code yet...")


## Make the output directories
def makeDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
def makeOutFileName(fastq):
	outdir = "/".join(fastq.split("/")[:-2]) + "/stats/"
	makeDir(outdir)
	outname = "raw_read_stats.dat"
	return outdir + outname

## Define some parameters
qualbins = range(60)

stats = {"Total2dBases":0}
header = ["readType","id","length","meanQualScore","sdQualScore","GC"]
logging.info("Analysing fastq files")
# outfile_name = makeOutFileName(fqFilelist[0])
outfile_name = opt.outfile_name
print outfile_name
with open(outfile_name,'w') as outfile:
	writer = csv.writer(outfile,delimiter="\t")
	writer.writerow(header)
	for fastq in fqFilelist:

		read_type = get_readfastq.split("/")[-1].split(".fq")[0]
		print readType
		## Storing values for histograms
		runningQualityScoreHist = [0]*len(qualbins)		
		## Write stats to file
		print fastq
		for record in SeqIO.parse(fastq, "fastq-sanger"):
			if readType == "2d":
				stats["Total2dBases"] += len(record.seq)
			qualityScores = record.letter_annotations["phred_quality"]
			row = [readType,record.id,len(record.seq),np.mean(qualityScores),np.std(qualityScores),GC(record.seq)]
			hist,bin_edges=np.histogram(qualityScores,bins=qualbins)
			runningQualityScoreHist = [i+j for i,j in zip(hist,runningQualityScoreHist)]
			writer.writerow(row)
print stats
for fast5 in fast5Filelist:
	f5 =  fast5File(filepath=f)
	


