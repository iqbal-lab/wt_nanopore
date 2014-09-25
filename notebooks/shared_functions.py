# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np

# <codecell>

def tabulate(x):
    u, indices = np.unique(x, return_inverse=True)
    return dict(zip(u, np.bincount(indices)))

