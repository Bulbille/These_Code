#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
dH1 = plt.subplot(111)
import numpy as np

L = 8

## Hamiltonien normal
for i in np.arange(L):
    for j in np.arange(L) :
        if j < L-1 :
            dH1.scatter(i,j,marker='o',color="blue")
            if i < L-1  :
                dH1.plot([i,i+1],[j,j],color="blue")
        if j < L-2 :
            dH1.plot([i,i],[j,j+1],color="blue")
    dH1.scatter(i,L,marker='o',color="blue")
dH1.plot([0,L-1],[L,L],color="blue")

dH1.arrow(-0.5,-0.5,L,0,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH1.arrow(L,-0.5,-L,0,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH1.text(L/2-0.5,-1.5,'$L$',color="blue",fontsize=20)

dH1.arrow(L,-0.5,0,L-1.5,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH1.arrow(L,L-1.5,0,-L+1.5,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH1.text(L+0.5,L/2-1.5,'$L\'-1$',color="blue",fontsize=20)
dH1.text(L/2-1,L-1.5,'$+$',color="blue",fontsize=20)

dH1.axis('off')
#dH1.axis('square')
dH1.set_aspect('equal')
dH1.set_xlim([-1,L+1])
dH1.set_ylim([-1,L+1])
plt.savefig('cross-h1.pdf')
plt.show()
