# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from shared_imports import *
from shared_functions import *

# <codecell>

all_reads = list(pysam.Samfile("/data/minion/work_old/expt05_aliquot_1_workflow_1_9_1/last/3D7_V3/reads_2D.last.sorted.bam", "rb" ).fetch())

# <codecell>

all_reads[0].aligned_pairs

# <codecell>


