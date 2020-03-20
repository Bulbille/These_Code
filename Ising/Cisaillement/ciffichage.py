#!/usr/bin/python
# -*- coding: utf8


import numpy as np
import matplotlib.pyplot as plt
fontsize= 12
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')
plt.rc("text",usetex=True) #Latex
import sys
import os
import re

#####################################
# Affiche l'état d'un système de 
# de taille N*M 
#####################################

S=np.loadtxt(sys.argv[1])
shear = re.findall("(\d\.\d+)",sys.argv[1])[0][:4]
M=int(S[-1][1])+1 # direction y 
N=int(S[-1][0])+1 #direction x 
s=np.zeros((M,N))

for k in range(len(S)):
    s[int(S[k,1]),int(S[k,0])]=int(S[k,2])

fig = plt.subplot(111)
fig.pcolormesh(s,rasterized=True,cmap="Pastel1")

print(shear,float(shear))

for y in np.arange(M):
    if (y+2) % 5 == 0:
        fig.arrow(N/2,y,(y*2/(M+1)-1)*30*float(shear),0,color="blue",length_includes_head=True,head_starts_at_zero=True,head_width=3,head_length=3)
fig.scatter(np.ones(M)*N/2,np.arange(M),s=1,color="blue")
fig.text(N/2+5,M/2-5,"\\textbf{$\omega="+shear+"$}",fontsize=20,color="blue")

fig.axis("image")
fig.set_ylim([0,M])
fig.set_xlabel('$x$')
fig.set_ylabel('$y$')

plt.savefig('cis-ising-f-'+shear+'.pdf', bbox_inches='tight')
plt.show()
plt.close()

