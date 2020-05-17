#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
dessin = plt.subplot(111)
import numpy as np
fs = 20
L = 8

#Définition de la ligne
x = np.ceil(np.arange(2*L)/2.)
a = [1,-1,-2,0,3,2,1,-2]
b = np.empty(2*len(a))
for j,i in enumerate(b):
    print(j,j/2,i)
    if j%2 == 0:
        b[j] = a[int(j/2)]+L/2
    else :
        b[j] = b[j-1]
print(a,b)
dessin.plot(x,b,color="black")

#Définition des + et -
pm = np.empty([L,L])
for i in np.arange(L) :
    plt.plot([i,i],[0,a[i]+L/2],color="black")
#    dessin.text(i+0.3,1,"-",fontsize=30)
plt.plot([L,L],[0,a[i]+L/2],color="black")
plt.fill_between(x,b,color="blue")
plt.plot([0,L],[0,0],color="black")
plt.plot(np.linspace(0,L,300),np.linspace(0,L,300)*0+L/2,'-.',color="red",linewidth=5)
#dessin.text(4.3,7.5,"+",fontsize=30)

#Légende des hauteurs
for i in np.arange(L) :
    dessin.text(i+0.5,-0.7,str(a[i]),fontsize=fs)

dessin.text(-1.3,-0.7,'$h_i = $',fontsize=fs)
dessin.text(-1.3,L/2-0.25,'$z =0 $',fontsize=fs)

#dessin.axes.get_xaxis().set_visible(False)
#dessin.axes.get_yaxis().set_visible(False)
dessin.axis('off')
dessin.set_ylim(-1,8)
dessin.set_xlim(-1,8)
plt.tight_layout()
plt.savefig('sos-indiscernable-inf.png')
plt.show()
