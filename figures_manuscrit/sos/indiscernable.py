#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
dessin = plt.subplot(111)
import numpy as np

L = 8

#Définition de la ligne
x = np.ceil(np.arange(2*L)/2.)
a = [4,5,3,4,7,6,4,2]
b = np.empty(2*len(a))
for j,i in enumerate(b):
    print(j,j/2,i)
    if j%2 == 0:
        b[j] = a[int(j/2)]
    else :
        b[j] = b[j-1]
dessin.plot(x,b)

#Définition des + et -
pm = np.empty([L,L])
for i in np.arange(L) :
    plt.plot([i,i],[0,a[i]],color="blue")
    dessin.text(i+0.3,1,"-",fontsize=30)
plt.plot([L,L],[0,a[i]],color="blue")
dessin.text(4.3,7.5,"+",fontsize=30)

#Légende des hauteurs
for i in np.arange(L) :
    dessin.text(i+0.5,-0.5,str(a[i]))

dessin.text(-0.5,-0.5,'$h_i = $',fontsize=16)

#dessin.axes.get_xaxis().set_visible(False)
#dessin.axes.get_yaxis().set_visible(False)
dessin.axis('off')
dessin.set_ylim(-1,8)
dessin.set_xlim(-1,8)
plt.savefig('sos-indiscernable.pdf')
plt.show()
