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
regex_chain="T-([\d.]+)(?:-\d)?(?:-h-([\d.]+))?(?:-time-([\d.]+))?"
directory = sys.argv[1]

for file in os.listdir(directory):
    regex = re.search(regex_chain,file)
    if regex == None :
        continue
    res = re.findall(regex_chain,file)
    time = res[0][2]
    if time == '':
        continue
    print file

    S=np.loadtxt(directory+file)
    M=int(S[-1][1])+1
    N=int(S[-1][0])+1
    s=np.zeros((M,N))
    for k in range(len(S)):
        s[int(S[k,1]),int(S[k,0])]=int(S[k,2])

    plt.figure()
    plt.pcolormesh(s,rasterized=True,cmap="Pastel1")
    plt.axis("tight")
    plt.title('Time = '+str(time))
    plt.savefig(directory+'temps-'+'{:0>7}'.format(time)+'.pdf')

import subprocess
command = "convert -delay 1x4 -loop 0 "+directory+"temps*.pdf "+directory+"animation.gif"
os.system(command)
