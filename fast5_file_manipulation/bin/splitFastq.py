#! /usr/bin/env python

## Takes a fastq files and splits each read into reads of length N

from optparse import OptionParser
from Bio import SeqIO

parser = OptionParser()
parser.description = """
Takes a fastq files and splits each read into reads of length N"""
parser.add_option("-f", "--file", dest="filepath",
                  help="path to fastq file")
parser.add_option("-n", "--N", dest="N",
                  help="length of outputted reads",default=1000,type="int")
(opt, args) = parser.parse_args()

output_handle = open(opt.filepath+".split%i" % opt.N, "w") 
sequences = []
for read in SeqIO.parse(opt.filepath,'fastq'):
	if len(read.seq) > opt.N:
		ran = range(1,len(read.seq),100)
		for i in range(len(ran)):
			tmpseq = read[ran[i]:ran[i]+opt.N]
			tmpseq.id += "_i_%i_len_%i"% (i,len(tmpseq))
			sequences.append(tmpseq)
	else:
		sequences.append(read)
SeqIO.write(sequences, output_handle, "fasta")