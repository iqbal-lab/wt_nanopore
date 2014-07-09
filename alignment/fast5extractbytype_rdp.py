#!/usr/bin/env python

# ============================================================================ #
#                                                                              #
# fast5extractbytype.py                                                        #
#                                                                              #
# ============================================================================ #

import h5py
from Bio import SeqIO
from StringIO import StringIO
import os, sys

if len(sys.argv) < 5:
    print 'Usage: {0} indir outreadtype outformat outdir'.format(os.path.basename(sys.argv[0]))
    print '       where outformat must be "fasta" or "fastq".'
    print ''
    sys.exit(1)
indir = os.path.expandvars(sys.argv[1])
outreadtype = sys.argv[2]
outformat = sys.argv[3]
outdir = os.path.expandvars(sys.argv[4])

log_path = os.path.join(outdir, '{0}_{1}.log'.format('.'.join(os.path.basename(sys.argv[0]).split('.')[0:-1]), outformat))
log_fp = open(log_path, 'w')

out_path = os.path.join(outdir, 'reads_{1}.{2}'.format('.'.join(os.path.basename(sys.argv[0]).split('.')[0:-1]), outreadtype, outformat))
out_fp = open(out_path, 'w')

keys = {'template' : '/Analyses/Basecall_2D_000/BaseCalled_template/Fastq',
        'complement' : '/Analyses/Basecall_2D_000/BaseCalled_complement/Fastq',
        '2D' : '/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'}

for subdir, dirs, files in os.walk(indir):
    for file in files:
        if file.endswith('.fast5'):
            fast5_path = os.path.join(subdir, file)
            hdf = h5py.File(fast5_path, 'r')
            numfound = 0
            for id, key in keys.iteritems():
                if id == outreadtype:
                    try:
                        fq = hdf[key][()]
                        run_id =  hdf['UniqueGlobalKey/tracking_id'].attrs['run_id']
                        numfound += 1
                        rec = SeqIO.read(StringIO(fq), 'fastq')
                        rec.id += "_" + id
                        rec.id = "run_id_" + run_id + '_' + rec.id
                        rec.description = 'len={0}'.format(len(rec.seq))
                        mode = 'w' if numfound == 1 else 'a'
                        SeqIO.write([rec], out_fp, outformat)
                    except Exception, e:
                        pass
                    if numfound == 0:
                        log_fp.write('Warn: File did not contain any template, complement or 2D reads ({0})\n'.format(fast5_path))
            hdf.close()

out_fp.close()
log_fp.close()

