#!/usr/bin/env python
import sys
sys.path.append('../src')
from fast5 import fast5File
import os
import glob
import logging
from optparse import OptionParser
parser = OptionParser()
parser.description = """
Converts an ONT fast5 file into fastq files saved to outdir/{template,complement,2d}"""
parser.add_option("-d", "--dir", dest="directory",
                  help="path to ONT data/reads/downloads directory")
parser.add_option("-o", "--out", dest="outdir",
                  help="path to output directory",default="./")
(opt, args) = parser.parse_args()
if opt.directory is None:
    parser.print_help()
## Make the output directories
def makeDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
makeDir(opt.outdir)

directory = opt.directory + "/*.fast5"
fast5list =  glob.glob(directory)
with open('complement.fq', 'w') as complementOutFile, \
     open('template.fq', 'w') as templateOutFile, \
     open('2d.fq', 'w') as twoDOutFile:
    for f in fast5list:
        fq =  fast5File(filepath=f)
        if fq.complement is not None:
            complementOutFile.write(fq.complement)
        if fq.template is not None:
            templateOutFile.write(fq.template)
        if fq.twoD is not None:
            twoDOutFile.write(fq.twoD)




