#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
from matplotlib import patches
dHlambda = plt.subplot(111)
import numpy as np

L = 8

## Hamiltonien normal
for i in np.arange(L):
    for j in np.arange(L) :
        dHlambda.scatter(i,j,marker='o',color="blue")
        if i < L-1 :
            dHlambda.plot([i,i+1],[j,j],color="blue")
        if j < L-1 :
            dHlambda.plot([i,i],[j,j+1],color="blue")



#### Rajout par dessus des liens destructeurs d'abord
k =4
dHlambda.text(L-0.2,k,'$k$',color="blue",       fontsize=12)
dHlambda.text(L-0.2,k-1,'$k-1$',color="blue",   fontsize=12)
dHlambda.text(L-0.2,k+1,'$k+1$',color="blue",   fontsize=12)

for i in np.arange(L):
    dHlambda.plot([i,i],[k,k+1],color="white")
    dHlambda.plot([i,i],[k,k-1],color="white")
    dHlambda.plot([i,i],[k,k+1],'-.',color="green")
    dHlambda.plot([i,i],[k,k-1],'-.',color="green")
    arc = patches.Arc((i,k),1,2,0,theta1=90,theta2=270,color="red",linestyle="--")
    dHlambda.add_patch(arc)

dHlambda.axis('off')
#dHlambda.axis('tight')
#dHlambda.set_aspect('equal')
#dHlambda.set_xlim([-1,L+1])
#dHlambda.set_ylim([-1,L+1])
plt.tight_layout(pad=-1.0)
plt.savefig('cross-hlambda.pdf',bbox_inches='tight')
plt.show()
