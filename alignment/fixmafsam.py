#!/usr/bin/env python

# ============================================================================ #
# fixmafsamfile.py                                                             #
# ============================================================================ #
'''
Fixes SAM file produced by lastal and maf-convert.py

The lastal program from the LAST package aligns reads to a reference and outputs
the alignment in MAF format. This can be converted to SAM format with the
maf-convert.py script that comes with the LAST package, but there are several
things wrong with the output.

This script outputs a SAM file with the following fixes:

1. Reads have base qualities in the correct strand orientation.
2. The file has a valid SAM header with a single read group.
3. Reads are placed into a single read group.
4. Hard-clipped (H) bases in the CIGAR string changed to soft-clipped (S).
   so CIGAR length (M+I+S) matches seq length in a valid SAM file.
4. Unmapped reads are appended to the end of the file and have:
   - mapping position of zero
   - mapping quality of zero
   - bitwise flag of 4
'''
# ============================================================================ #
# Camilla Ip                                                                   #
# camilla.ip@well.ox.ac.uk                                                     #
# Wellcome Trust Centre for Human Genetics, University of Oxford, UK           #
# June 2014                                                                    #
# ============================================================================ #

# ============================================================================ #
# Import Modules                                                               #
# ============================================================================ #

import argparse, os, shlex, subprocess as sp, sys
from subprocess import PIPE
from Bio import SeqIO

# ============================================================================ #
# Global variables                                                             #
# ============================================================================ #

_progdesc = 'Fixes SAM file produced by lastal and maf-convert.py'

# ============================================================================ #
# General-purpose functions                                                    #
# ============================================================================ #

def sys_exec(self, cmd):
    """
    Execute a command using the subprocess module to trap the
    return value, the stdout string and stderr string.
    """
    proc_handle = sp.Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
    proc_stdout, proc_stderr = proc_handle.communicate()
    proc_returncode = proc_handle.returncode
    return [proc_returncode, proc_stdout, proc_stderr]

# ============================================================================ #
# Program usage                                                                #
# ============================================================================ #

def Initialise():
    'Read in the command-line arguments. Check for errors.'

    # Process command-line arguments
    global _progdir, _progname, _args
    _progdir = os.path.dirname(os.path.realpath(sys.argv[0]))
    _progname = os.path.basename(os.path.realpath(sys.argv[0]))
    progexamples = []
    parser = argparse.ArgumentParser(description=_progdesc, \
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fastqreads', metavar='str', dest='fastqreads', \
        type=str, default=None, help='Path to reads in fastq format', required=True)
    parser.add_argument('--lastmaf', metavar='str', dest='lastmaf', \
        type=str, default=None, help='Path to alignments from lastal in MAF format', required=True)
    parser.add_argument('--lastsam', metavar='str', dest='lastsam', \
        type=str, default=None, help='Path to reads from lastal and maf-convert.py in SAM in format', required=True)
    parser.add_argument('--readgroup', metavar='str', dest='readgroup', \
        type=str, default=None, help='The read group that all reads will be assigned to', required=True)
    parser.add_argument('--fastaref', metavar='str', dest='fastaref', \
        type=str, default=None, help='Path to reference contig in fasta format ' \
        '(only works for references with 1 contig)', required=True)
    parser.add_argument('--outsampath', metavar='str', dest='outsampath', \
        type=str, default=None, help='Path for output SAM file', required=True)
    _args = parser.parse_args()

    # Error checking
    if not os.path.exists(os.path.expandvars(_args.fastqreads)):
        sys.stderr.write('Error: fastq reads file does not exist ({0})\n'.format(_args.fastqreads))
        sys.exit(2)
    if not os.path.exists(os.path.expandvars(_args.lastsam)):
        sys.stderr.write('Error: maf-sam reads file does not exist ({0})\n'.format(_args.lastsam))
        sys.exit(2)
    if os.path.exists(os.path.expandvars(_args.outsampath)):
        sys.stdout.write('Warn: outsampath already exists - overwriting ({0})\n'.format(_args.outsampath))

# ============================================================================ #
# Do the fixing                                                                #
# ============================================================================ #

def Get_ReferenceLength(refpath):
    'Return the name and length of the reference genome'
    data = SeqIO.read(os.path.expandvars(refpath), 'fasta')
    return str(data.id), len(data.seq)

def Get_SeqAndBaseQualitiesForEachRead(fastqpath):
    'Parse the fastq file of reads, save the bq string in a dictionary keyed on readid.'
    D = {}
    in_fp = open(fastqpath, 'r')
    if not in_fp:
        sys.stdout.stderr('Error: failed to open fastqreads file ({0})\n'.format(_args.fastqreads))
        sys.exit(3)
    for record in SeqIO.parse(in_fp, 'fastq'):
        #print record.id, record.seq[0:10], ''.join([chr(x+33) for x in record.letter_annotations.values()[0][0:10]])
        if not D.has_key(record.id):
            D[record.id] = [str(record.seq), ''.join([chr(x+33) for x in record.letter_annotations.values()[0]])]
        else:
            sys.stderr.write('Warn: more than one record with same id - overwriting previous bq value ({0})\n'.format(record.id))
    in_fp.close()
    return D

def  Get_BestLastAlnForEachRead(lastmafpath):
    '''
    Parse the lastal output file in MAF format and save the best mapping score
    and the 1-based position at which the alignment starts for that readid.
    '''
    bestscore = {}
    aln = { 'score':0, 'name':0, 'start':0, 'alnsize':0, 'strand':0, 'seqsize':0, 'alignment':0 }
    with open(os.path.expandvars(lastmafpath), 'r') as in_fp:
        for line in in_fp:
            line = line.strip()
            if line.startswith('#'):
                continue
            if line.startswith('a'):	# Have score line
                blockcnt = 1
                aln = dict.fromkeys(aln, 0)
                aln['score'] = int(line.split('=')[1])
                continue
            if blockcnt == 1:
                blockcnt += 1
                aln['start'] = int(line.split()[2]) + 1
                continue
            if blockcnt == 2:
                blockcnt += 1
                L = line.split()
                #aln['name'], aln['start'], aln['alnsize'], aln['strand'], aln['seqsize'], aln['alignment'] = L[1:]
                aln['name'] = L[1]
                aln['alnsize'] = L[3]
                aln['strand'] = L[4]
                #aln['start'] = int(aln['start']) + 1
                aln['alnsize'] = int(aln['alnsize'])
                if not bestscore.has_key(aln['name']):
                    bestscore[aln['name']] = [int(aln['score']), int(aln['start']), int(aln['alnsize'])]
                else:
                    if bestscore[aln['name']][0] < aln['score']:
                        if bestscore[aln['name']][2] > aln['alnsize']:
                            sys.stdout.write('Warn: replacing bestscore with shorter alignment (readid={0}: score={1}->{2}, alnsize={3}->{4})\n'.format( \
                                aln['name'], \
                                bestscore[aln['name']][0], aln['score'], \
                                bestscore[aln['name']][2], aln['alnsize']))
                        bestscore[aln['name']] = [int(aln['score']), int(aln['start']), int(aln['alnsize'])]
                blockcnt = 0
    return bestscore

def Print_NewSAMFile(seqandbq, bestscore, refname, reflen, mafsampath, fastqpath, readgroup, outsampath):
    '''
    Read each line of mafsampath and output the corrected SAM records in the
    same order. Then append the unmapped reads at the end in any order.
    Since it is possible for there to be more than one alignment with the same
    score and length, following preliminary spot-check to verify that neither
    alignment is very different (usually just in exact position of a gap in
    regions that are difficult to align), a dictionary of readids printed is
    retained and only the first is printed to the SAM file.
    '''
    with open(outsampath, 'w') as out_fp:
        readidmapped = []
        # Fix all the mapped reads and create a list of all the readids for mapped reads
        with open(mafsampath, 'r') as in_fp:
            recordcnt = 0
            for line in in_fp:
                if line.startswith('@'):
                    # Skip header lines
                    out_fp.write('{0}'.format(line))
                else:
                    # Add additional header lines
                    recordcnt += 1
                    oldL = line.strip().split('\t')
                    if recordcnt == 1:
                        out_fp.write('{0}\n'.format('\t'.join(['@SQ', 'SN:{0}'.format(refname), 'LN:{0}'.format(reflen)])))
                        out_fp.write('{0}\n'.format('\t'.join(['@RG', 'ID:{0}'.format(readgroup)])))
                    readid = oldL[0]
                    if readid == 'channel_386_read_63_2D':
                        pass
                    linerec_startpos = int(oldL[3])
                    bestrec_startpos = bestscore[readid][1]
                    if bestrec_startpos == linerec_startpos:
                        if readid not in readidmapped:
                            # Print if this is the best match for this readid
                            if seqandbq.has_key(readid):
                                if bestscore.has_key(readid):
                                    newmapqual = bestscore[readid][0]
                                else:
                                    newmapqual = oldL[5]
                                    sys.stdout.write('Warn: failed to retrieve last mapping score so using default supplied by maf-convert.py ({0}: mapqual={1})\n'.format(readid), newmapqual)
                                #newL = oldL[0:4] + [newmapqual] + oldL[5:9] + seqandbq[readid] + oldL[11:] + ['RG:Z:{0}'.format(readgroup)]
                                newL = oldL[:]
                                newL[5] = newL[5].replace('H', 'S')
                                newL[4] = newmapqual
                                newL[9] = seqandbq[readid][0]
                                newL[10] = seqandbq[readid][1]
                                newL += ['RG:Z:{0}'.format(readgroup)]
                                # Debugging
                                #newL[9] = newL[9][0:10]
                                #newL[10] = newL[10][0:10]
                                out_fp.write('{0}\n'.format('\t'.join([str(x) for x in newL])))
                            else:
                                sys.stdout.write('Warn: failed to find base qualities for read with id ({0})\n'.format(readid))
                            if readid in readidmapped:
                                sys.stdout.write('Warn: adding readid to list of mapped ids again ({0})\n'.format(readid))
                            readidmapped.append(readid)
                        else:
                             sys.stdout.write('Warn: skipping line because already kept the 1st match with same starting position and score ({0}: start={1})\n'.format(readid, oldL[3]))
                    else:
                        sys.stdout.write('Warn: skipping line because it is not the best LAST match ({0}: start={1})\n'.format(readid, oldL[3]))
        # Output all the unmapped reads in any order.
        readidunmapped = set(seqandbq.keys()).difference(set(readidmapped))
        for readid in set(readidmapped):
            print 'Debug: Mapped: {0}'.format(readid)
        sys.stdout.write('Info: {0} mapped reads, {1} unmapped reads\n'.format(len(readidmapped), len(readidunmapped)))
        in_fp = open(fastqpath, 'r')
        if not in_fp:
            sys.stdout.stderr('Error: failed to open fastqreads file ({0})\n'.format(_args.fastqreads))
            sys.exit(3)
        for record in SeqIO.parse(in_fp, 'fastq'):
            if record.id in readidunmapped:
               newL = [record.id, 4, refname, 0, 0, '*', '=', 0, 0, \
                   str(record.seq), seqandbq[record.id][1], \
                   'MQ:i:0', 'RG:Z:{0}'.format(readgroup)]
               out_fp.write('{0}\n'.format('\t'.join([str(x) for x in newL])))
        in_fp.close()

def Fix_SAM_File():
    ''
    refname, reflen = Get_ReferenceLength(_args.fastaref)
    seqandbq = Get_SeqAndBaseQualitiesForEachRead(_args.fastqreads)
    bestscore = Get_BestLastAlnForEachRead(_args.lastmaf)
    Print_NewSAMFile( \
        seqandbq, \
        bestscore, \
        refname, reflen, \
        os.path.expandvars(_args.lastsam), \
        os.path.expandvars(_args.fastqreads), \
        os.path.expandvars(_args.readgroup), \
        os.path.expandvars(_args.outsampath))

# ============================================================================ #
# Main                                                                         #
# ============================================================================ #

if __name__ == '__main__':

    Initialise()
    Fix_SAM_File()
    sys.exit(0)
    try:
        Initialise()
        Fix_SAM_File()
    except:
        sys.stderr.write('{0}\n'.format('\n'.join([str(x) for x in sys.exc_info()])))
        sys.exit(4)

# ============================================================================ #
