#!/usr/bin/python
# -*- coding: utf8


import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import re

#####################################
# Affiche l'état d'un système de 
# de taille N*M 
#####################################
fig, ax = plt.subplots(nrows=2, ncols=2)

ttc = "0.70"
temps = [1]
i = 0
for row in ax:
    for col in row:
        S=np.loadtxt('glau-interface/snap'+ttc+'_time'+str(temps[i]))
        M=int(S[-1][1])+1
        N=int(S[-1][0])+1
        s=np.zeros((M,N))

        for k in range(len(S)):
            s[int(S[k,1]),int(S[k,0])]=int(S[k,2])
        col.pcolormesh(s,rasterized=True,cmap="Pastel1")
        col.set_title("t="+str(10*temps[i])+' MC')
        col.axis('off')
        i+= 1

#plt.figure()
#plt.axis("image")
#plt.xlabel("X")
#plt.ylabel("Y")

#try :
#	interface=np.loadtxt(sys.argv[1]+"-interface")
#	plt.plot(interface,'o')
#except :
#	pass

plt.savefig('clusterization.pdf')
plt.show()
plt.close()
