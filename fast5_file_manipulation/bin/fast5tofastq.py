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
Converts an ONT fast5 file into fastq files saved to outdir/{template,complement,2d}.fq"""
parser.add_option("-d", "--dir", dest="directory",
                  help="path to ONT data/reads/downloads directory")
parser.add_option("-o", "--out", dest="outdir",
                  help="path to output directory",default="./")
(opt, args) = parser.parse_args()
if opt.directory is None:
    parser.print_help()

def create_dir(path):
  if not os.path.exists(path):
      os.makedirs(path)  

create_dir(opt.outdir)
create_dir('%s/complement' % opt.outdir)  
create_dir('%s/template' % opt.outdir)  
create_dir('%s/2d' % opt.outdir)  


directory = opt.directory + "/*.fast5"
fast5list =  glob.glob(directory)
# with open('%s/complement.fq' % opt.outdir, 'w') as complementOutFile, \
#      open('%s/template.fq' % opt.outdir, 'w') as templateOutFile, \
#      open('%s/2d.fq' % opt.outdir, 'w') as twoDOutFile:
#     for f in fast5list:
#         try:

#           fq =  fast5File(filepath=f)
#           if fq.complement is not None:
#               complementOutFile.write(fq.complement)
#           if fq.template is not None:
#               templateOutFile.write(fq.template)
#           if fq.twoD is not None:
#               twoDOutFile.write(fq.twoD)
#           fq.close()
#         except IOError:
#           print "Passing %s " % f 


for f in fast5list:
  basename =  os.path.basename(f).split('.')[0]
  try:
    fq =  fast5File(filepath=f)
    if fq.complement is not None:
      with open('%s/complement/%s_complement.fq' % (opt.outdir,basename), 'w') as complementOutFile:
        complementOutFile.write(fq.complement)
    if fq.template is not None:
      with open('%s/template/%s_template.fq' % (opt.outdir,basename), 'w') as templateOutFile:      
        templateOutFile.write(fq.template)
    if fq.twoD is not None:
      with open('%s/2d/%s_2d.fq' % (opt.outdir,basename), 'w') as twoDOutFile:      
        twoDOutFile.write(fq.twoD)   
    fq.close()  
  except (IOError,ValueError), e:
     print "Passing %s : %s " % (f ,e)  








