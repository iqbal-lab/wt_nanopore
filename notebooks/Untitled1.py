# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from shared_imports import *
from shared_functions import *

# <codecell>

temp = np.array(["a", "2"])

# <codecell>

''.join(["a", "2"])

# <codecell>

temp

# <codecell>

''.join(temp)

# <codecell>

import time
def convert_time(YOUR_EPOCH_TIME):
    return time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(YOUR_EPOCH_TIME))
convert_time(1403178270)

# <codecell>

'___'.join(map(convert_time, np.array([1403178270, 1403170000])))

# <codecell>


