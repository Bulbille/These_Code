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
    for j in np.arange(L) :
        if j >= a[i] :
            dessin.scatter(i+0.5,j+0.5,marker='+',s=80,color="blue")
        else :
            dessin.scatter(i+0.5,j+0.5,marker='_',s=80,color="blue")

#Légende des hauteurs
for i in np.arange(L) :
    dessin.text(i+0.5,-0.5,str(a[i]))
dessin.text(-0.5,-0.5,'$h_i = $',fontsize=16)

#dessin.axes.get_xaxis().set_visible(False)
#dessin.axes.get_yaxis().set_visible(False)
dessin.axis('off')
dessin.set_ylim(-1,8)
dessin.set_xlim(-1,8)
plt.savefig('figure-sos.pdf')
plt.show()
