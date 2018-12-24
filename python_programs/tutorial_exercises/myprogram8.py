#!usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

filename = '/Users/md875/data/tutorial_data/S1998148172338_chlor_a.f999x999'
f = open(filename)
full_image = np.fromfile(f, dtype = np.float32)
f.close()
full_image = full_image.reshape([999,999])

sub_image = full_image[300:700, 100:400]

#imaging
myColorScale = plt.get_cmap('spectral')
myColorScale.set_bad('k')

plt.figure(1)
plt.imshow(full_image, cmap = myColorScale, vmin = 0.01, vmax = 20.0, norm = matplotlib.colors.LogNorm())
plt.figure(2)
plt.imshow(sub_image, cmap = myColorScale, vmin = 0.01, vmax = 20.0, norm = matplotlib.colors.LogNorm())

plt.show()