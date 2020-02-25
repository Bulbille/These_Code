#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fontsize=14
plt.rc("font",size=fontsize)

dessin = plt.subplot(111)
import numpy as np

L = 6

a = ['+','+','-','-','+','+']
b = ['+','+','+','-','+','+']
c = ['+','-','-','-','+','+']
#Légende des hauteurs
for i in np.arange(L) :
    dessin.text(i,0,str(a[i]))
    dessin.text(i,0.5,str(b[i]))
    dessin.text(i,1,str(c[i]))
dessin.text(L,0,str('$M=2$'),color='blue')
dessin.text(L,0.5,str('$M=4$'),color='red')
dessin.text(L,1,str('$M=0$'),color='red')

#dessin.axes.get_xaxis().set_visible(False)
#dessin.axes.get_yaxis().set_visible(False)
dessin.axis('off')
dessin.set_ylim(-0,1)
dessin.set_xlim(-1,6)
plt.savefig('figure-canonique.pdf')
plt.show()
