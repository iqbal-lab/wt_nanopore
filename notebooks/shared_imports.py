# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
import vcfnp
import time
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from StringIO import StringIO

import vcf
from vcf.parser import _Info as Info
import vcfnp
print 'vcfnp', vcfnp.__version__
import vcfplt
print 'vcfplt', vcfplt.__version__


#import vcfarray
import petl.interactive as etl
from petl.interactive import *
import scipy
from scipy import stats
import numpy as np
from matplotlib import *
from matplotlib.pylab import *
import matplotlib.pyplot as plt
from numpy.lib.recfunctions import append_fields
import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
import pysam
import numexpr
import pysamstats
import petl
print 'petl', petl.VERSION
import petlx
print 'petlx', petlx.VERSION
import petlx.vcf
import petlx.interval
import petlx.gff3
import petlx.xlsx
from petlx.xlsx import *
from petlx.xls import *
from petlx.array import *
import petlx.ipython
import h5py
import tables
from datetime import datetime
from dateutil.relativedelta import relativedelta
from petlx.tabix import fromtabix
from bisect import bisect_left, bisect_right
import scipy.cluster.hierarchy as sch
from collections import OrderedDict
from IPython.display import display
import brewer2mpl
import subprocess

print 'numpy', np.__version__
print 'numexpr', numexpr.__version__
print 'pysam', pysam.__version__
print 'pysamstats', pysamstats.__version__
print 'petl', petl.VERSION
print 'petlx', petlx.VERSION
print 'vcf', vcf.VERSION
print 'vcfnp', vcfnp.__version__
print 'h5py', h5py.__version__
print 'tables', tables.__version__

