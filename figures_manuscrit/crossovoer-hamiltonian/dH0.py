#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
dH0 = plt.subplot(111)
import numpy as np

L = 8

## Hamiltonien normal
for i in np.arange(L):
    for j in np.arange(L) :
        dH0.scatter(i,j,marker='o',color="blue")
        if i < L-1 :
            dH0.plot([i,i+1],[j,j],color="blue")
        if j < L-1 :
            dH0.plot([i,i],[j,j+1],color="blue")

dH0.arrow(-0.5,-0.5,L,0,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH0.arrow(L,-0.5,-L,0,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH0.text(L/2-0.5,-1.5,'$L\'$',color="blue",fontsize=20)

dH0.arrow(L,-0.5,0,L-0.5,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH0.arrow(L,L,0,-L,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH0.text(L+0.5,L/2-0.5,'$L$',color="blue",fontsize=20)

dH0.axis('off')
#dH0.axis('square')
#dH0.set_aspect('equal')
plt.tight_layout(pad=-1.0) 
dH0.set_xlim([-1,L+1])
dH0.set_ylim([-1,L+1])
plt.savefig('cross-h0.pdf',bbox_inches='tight')
plt.show()
