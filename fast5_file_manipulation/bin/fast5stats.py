#!/usr/bin/env python
import sys
sys.path.append('../src')
from fast5 import fast5File
import glob
import logging
import csv
from optparse import OptionParser

parser = OptionParser()
parser.description = """
Looks at fast5 files in directory and returns stats on comtents of files"""
parser.add_option("-d", "--dir", dest="directory",
                  help="path to ONT data/reads/downloads directory")
parser.add_option("-o", "--out", dest="outfilepath",
                  help="path to output file",default=None)
(opt, args) = parser.parse_args()

if opt.directory is None:
    parser.print_help()

def printStats(stats,last=False):
	template_out = "%i template " % stats.get('templateReads',0) 
	complement_out = "%i complement" % stats.get('complementReads',0) 
	twod_out = "%i 2d reads out of " % stats.get('2dReads',0) 
	# noFq_out = "%i of fast5 files without fq read\r" % stats.get('noFq',0) 
	total_out = "%i fast5 files\r" % stats.get('total',0) 
	if not last:
		sys.stdout.write(" ".join([template_out,complement_out,twod_out,total_out]))
		sys.stdout.flush()
	else:
		sys.stdout.write("\n".join(["\n"+template_out,complement_out,twod_out,total_out+'\n']))
	
	# print "%i of fast5 files with complement reads without template reads" % stats.get('complementReadsWithoutTemplate',0)

directory = opt.directory + "/*.fast5"
fast5list =  glob.glob(directory)

stats = {'templateReads':0,
		'complementReads':0,
		'2dReads':0,
		'complementReadsWithoutTemplate':0,
		'noFq':0,
		'total':0}
if opt.outfilepath:
	outfile = open(opt.outfilepath, 'w')
	writer = csv.writer(outfile,delimiter="\t")
	header = ["filename","has2dfq","hasComplementfq","hasTemplatefq"]
	writer.writerow(header)

for f in fast5list:
    fq =  fast5File(filepath=f)
    has2d = int(fq.has2d)
    hasComplement = int(fq.hasComplement)
    hasTemplate = int(fq.hasTemplate)
    fq.close()
    if opt.outfilepath:
    	writer.writerow([f.split('/')[-1].split('.')[0],has2d,hasComplement,hasTemplate])
    stats['templateReads'] += hasTemplate
    stats['complementReads'] += hasComplement
    stats['2dReads'] += has2d
    stats['complementReadsWithoutTemplate'] += (hasComplement and not hasTemplate)
    stats['noFq'] += (not hasComplement and not hasTemplate and not has2d)
    stats['total'] += 1
    printStats(stats)
if opt.outfilepath:
	outfile.close()
printStats(stats,last=True)	
	    

