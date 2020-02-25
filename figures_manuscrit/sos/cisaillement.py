#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fontsize=14
plt.rc("font",size=fontsize)

dessin = plt.subplot(111)
import numpy as np

LY = 10
LX = 10

x   = np.linspace(0,LX,200)
interface = 2*np.cos(x)+LY/2
dessin.plot(x,interface)

phi1 = 0.6
phi2 = 0.2

#Légende des hauteurs
for j in np.arange(LY) :
    for i in np.arange(LX) :
        if j < 2*np.cos(i)+LY/2 :
            dessin.arrow(i,j,phi1,0,color="red",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.3,head_length=0.1)
        else :
            dessin.arrow(i,j,phi2,0,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.3,head_length=0.1)
#    dessin.scatter(0,j)
#    dessin.scatter(phi1,j)

#dessin.axes.get_xaxis().set_visible(False)
#dessin.axes.get_yaxis().set_visible(False)
dessin.set_xlabel('$\\vec{x}$')
dessin.set_ylabel('$z$')
dessin.set_xlim(0,LX)
dessin.set_ylim(-0.2,LY-0.8)
plt.savefig('ab-phi.pdf')
plt.show()
