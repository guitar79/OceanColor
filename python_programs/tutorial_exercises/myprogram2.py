#!usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
yr = np.array(range(16)) + 2000
chl = np.random.rand(16)

plt.plot(yr,chl)
plt.ylabel ('Annual Chlorophyll'); plt.xlabel('Year');

plt.show()
plt.close()
