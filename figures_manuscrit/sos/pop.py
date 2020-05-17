#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
dessin = plt.subplot(111)
import numpy as np
fs = 15
L = 8

#Définition de la ligne
x = np.ceil(np.arange(2*L)/2.)

### Particles A
a = [4,5,3,4,7,6,4,2]
b = np.empty(2*len(a))
for j,i in enumerate(b):
    if j%2 == 0:
        b[j] = a[int(j/2)]
    else :
        b[j] = b[j-1]
dessin.plot(x,b,color="blue")

### Particles B
c = [1,2,4,4,0,3,2,3]
h = np.add(a,c)
d = np.empty(2*len(c))
for j,i in enumerate(b):
    if j%2 == 0:
        d[j] = c[int(j/2)]+a[int(j/2)]
    else :
        d[j] = d[j-1]
dessin.plot(x,d,color="green")
dessin.fill_between(x,d,color='lightgreen')
dessin.fill_between(x,b,color='mediumturquoise')

### Dessin des particules

for i in np.arange(L):
    for j in np.arange(max(h)):
        if j < a[i] :
            dessin.text(i+0.35,j+0.25,'$p_1$',fontsize=fs)
        elif j < h[i]:
            dessin.text(i+0.35,j+0.25,'$p_2$',fontsize=fs)

### Lignes de séparation verticles
for i in np.arange(L):
    dessin.plot([i,i],[0,a[i]],color="blue")
    dessin.plot([i,i],[a[i],h[i]],color="green")
dessin.plot([L,L],[0,a[-1]],color="blue")
dessin.plot([L,L],[a[-1],h[-1]],color="green")

#Légende des hauteurs
for i in np.arange(L) :
    dessin.text(i+0.35,-0.7,str(h[i]),fontsize=fs)
    dessin.text(i+0.35,-1.7,str(a[i]),fontsize=fs)
    dessin.text(i+0.35,-2.7,str(c[i]),fontsize=fs)

dessin.text(-1.2,-0.7,'$h_i = $',fontsize=fs)
dessin.text(-1.2,-1.7,'$n_{1,i} = $',fontsize=fs)
dessin.text(-1.2,-2.7,'$n_{2,i} = $',fontsize=fs)
dessin.text(-1.3,-0.1,'$z =0 $',fontsize=fs)

plt.plot(np.linspace(0,L,300),np.linspace(0,L,300)*0,'-.',color="red",linewidth=5)

#dessin.axes.get_xaxis().set_visible(False)
#dessin.axes.get_yaxis().set_visible(False)
dessin.axis('off')
dessin.set_ylim(-1,max(h))
dessin.set_xlim(-1,8)
plt.tight_layout()
plt.savefig('figure-pop.pdf')
#plt.show()
