#!/usr/bin/python
# -*- coding : utf8

import numpy as np
import matplotlib.pyplot as plt 
fontsize=12
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')
plt.rc("text",usetex=True) #Latex
import sys


data = np.loadtxt(sys.argv[1])

cercle1 = data[0:10]

from ellipses import LSqEllipse
from matplotlib.patches import Ellipse
lsqe = LSqEllipse()
lsqe.fit(cercle1)
center, width, height, phi = lsqe.parameters()

plt.close('all')
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.axis('equal')
ax.plot(data[0], data[1], 'ro', label='test data', zorder=1)

ellipse = Ellipse(xy=center, width=2*width, height=2*height, angle=np.rad2deg(phi),
                       edgecolor='b', fc='None', lw=2, label='Fit', zorder = 2)
ax.add_patch(ellipse)

plt.legend()
plt.show()
