#! /usr/bin/env python
## A script that takes a sam file as an argument and generates some basic stats about alignment

import sys
import os
sys.path.append('../src')
from optparse import OptionParser
import pysam
import logging
import csv
parser = OptionParser()
parser.description = """
Calculates some basic stats on ONT read alignment"""
parser.add_option("-f", "--file", dest="file",
                  help="path to sam/bam file")
# parser.add_option("-o", "--out", dest="outdir",
#                   help="path to output directory",default="./")
(opt, args) = parser.parse_args()
if opt.file is None:
    parser.print_help()
    sys.exit()
def makeDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
def makeOutFileName(infile,end):
    outdir = "/".join(infile.split("/")[:-2]) + "/stats/"
    makeDir(outdir)
    outname = infile.split("/")[-1].split(".")[0] +"_"+ end
    return outdir + outname


samfile = pysam.Samfile( opt.file )

covOutfile = makeOutFileName(opt.file,"coverage.dat")
print "Calculating coverage stats. Outputting to %s " % covOutfile
with open(covOutfile,'w') as outcov:
    writer = csv.writer(outcov,delimiter="\t")
    writer.writerow(["pos","cov"])
    for pileupcolumn in samfile.pileup():
        row = [pileupcolumn.pos , pileupcolumn.n]
        writer.writerow(row)


alignOutfile = makeOutFileName(opt.file,"alignment.dat")
print "Calculating alignment stats. Outputting to %s " % alignOutfile
with open(alignOutfile,'w') as outali:
    writer = csv.writer(outali,delimiter="\t")
    writer.writerow(["id","is_unmapped","pos","len","aligned_length","M","I","D","S","H"])
    for read in samfile.fetch():
        readLength = 0
        M,I,D,S,H = (0,0,0,0,0)
        for t,i in read.cigar:
            readLength += i
            if t == 0:
                ## Match or mismatch
                M += i 
            elif t == 1:                
                ## Insertion to the reference
                I += i
            elif t == 2:
                ## Deletion from the reference
                D += i 
            elif t == 4:
                ##  soft clipping (clipped sequences present in SEQ)
                S += i 
            elif t == 5:
                ##  hard clipping (clipped sequences NOT present in SEQ)
                H += i 
            else:
                pass
        row = [read.qname,read.is_unmapped,read.pos,readLength,read.alen,M,I,D,S,H]
        writer.writerow(row)


