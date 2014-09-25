#!/usr/bin/env python
#Konrad Paszkiewicz 2014, Exeter Sequencing Service
from StringIO import StringIO
import sys
import math
from collections import defaultdict

ids = {'Template only reads' : 'Splitter indicates no complement data. Workflow cannot continue',
        'Missing Template data' : 'Splitter indicates no template data. Workflow cannot continue',
        'Too few events' : 'Read has too few events. Workflow cannot continue',
        '2D basecalling successful' : 'Basecalling 2D hairpin data',
        'Poor ratio of template to complement data' : 'Ratio of template to complement data is no good. Workflow cannot continue',
        'Inconsistent template vs complement events' : 'Number of template events out of range aligning template',
        'Inconsistent complement vs template events' : 'Number of complement events out of range aligning'}

freqs={}
for i,v in ids.iteritems():
        freqs[i]=0


for filename in sys.argv[1:]:
        f = open(filename,'r')
        filecontents = f.read()
        flag=0
        for i,v in ids.iteritems():
                if v in filecontents:
                        freqs[i]=freqs[i]+1
                        flag=1
        if flag==0:
                print filename
total=0
for i,v in freqs.iteritems():
        total = total + v
for i,v in freqs.iteritems():
        percent = round(100*float(v)/float(total),2)
        print "%s  %i (%d percent)" % (i,v,percent)

print "Total    %i" % (total)