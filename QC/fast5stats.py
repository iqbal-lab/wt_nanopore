#!/usr/bin/env python

# ============================================================================ #
#                                                                              #
# fast5stats.py                                                                #
#                                                                              #
#                                                                              #
# ============================================================================ #

import h5py
from Bio import SeqIO
from StringIO import StringIO
import os, sys

if len(sys.argv) < 4:
    print 'Usage: {0} indir samplename outpath'.format(os.path.basename(sys.argv[0]))
    print 'Print statistics on all the fast5 files in indir.'
    print ''
    sys.exit(1)
indir = os.path.expandvars(sys.argv[1])
samplename = sys.argv[2]
outpath = sys.argv[3]

keys = {'template' : '/Analyses/Basecall_2D_000/BaseCalled_template/Fastq',
        'complement' : '/Analyses/Basecall_2D_000/BaseCalled_complement/Fastq',
        '2D' : '/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'}

H = ['samplename', 'channel', 'readnum', 'readtype', 'seqlen']
stats_fp = open(outpath, 'w')
stats_fp.write('{0}\n'.format('\t'.join(H)))

for subdir, dirs, files in os.walk(indir):
    for file in files:
        if file.endswith('.fast5'):
            fast5_path = os.path.join(subdir, file)
            hdf = h5py.File(fast5_path, 'r')
            readcount = 0
            for id, key in keys.iteritems():
                try:
                    fq = hdf[key][()]
                    readcount += 1
                    rec = SeqIO.read(StringIO(fq), 'fastq')
                    L = rec.id.split('_')
                    channelnum = L[1]
                    readnum = L[3]
                    readtype = id
                    seqlen = len(rec.seq)
                    result = [samplename, channelnum, readnum, readtype, seqlen]
                    stats_fp.write('{0}\n'.format('\t'.join([str(x) for x in result])))
                except Exception, e:
                    pass
            if readcount == 0:
                sys.stdout.write('Warn: No reads found ({0})\n'.format(fast5_path))
            hdf.close()

stats_fp.close()

