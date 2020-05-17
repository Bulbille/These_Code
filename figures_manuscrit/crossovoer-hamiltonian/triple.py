#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt 
fig,(dH0,dHlambda,dH1) = plt.subplots(nrows=1,ncols=3,figsize=(30,30))

from matplotlib import patches 
import numpy as np

L = 8 

## Hamiltonien 0
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


## Hamiltonien 1
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
dH1.text(L/2-0.5,-1.5,'$L\'$',color="blue",fontsize=20)

dH1.arrow(L,-0.5,0,L-1.5,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH1.arrow(L,L-1.5,0,-L+1.5,color="blue",shape="full",length_includes_head=True,head_starts_at_zero=True,head_width=0.5,head_length=1)
dH1.text(L+0.5,L/2-1.5,'$L-1$',color="blue",fontsize=20)
dH1.text(L/2-1,L-1.5,'$+$',color="blue",fontsize=20)


############### Hamitlonien crossover

for i in np.arange(L):
    for j in np.arange(L) :
        dHlambda.scatter(i,j,marker='o',color="blue")
        if i < L-1 :
            dHlambda.plot([i,i+1],[j,j],color="blue")
        if j < L-1 :
            dHlambda.plot([i,i],[j,j+1],color="blue")
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

dH0.axis('off')
dH0.set_xlim([-2,L+2])
dH0.set_ylim([-2,L+2])
dH0.axis('square')

dH1.axis('off')
dH1.set_xlim([-2,L+2])
dH1.set_ylim([-2,L+2])
dH1.axis('square')

dHlambda.axis('off')
dHlambda.axis('square')
dHlambda.set_xlim([-2,L+2])
dHlambda.set_ylim([-2,L+2])
#plt.tight_layout()

plt.savefig('hamilton-crossover.pdf',bbox_inches='tight')
plt.show()
