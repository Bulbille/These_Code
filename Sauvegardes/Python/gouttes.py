#!/user/bin/python
# -*- coding: utf8

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import scipy as sy
import numpy as np
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import fnmatch
import re

M=60
N=40


directory='./vrai/'


###Températures
temp = []
for file in sorted(os.listdir(directory)):
        if re.search("ttc_(\d+\.?\d*)",file) == None:
                continue
	t=re.findall("ttc_(\d+\.?\d*)",file)[0]
	if t not in temp:
		temp = np.append(temp,t)

###Plots
for t in temp:
	for file in sorted(os.listdir(directory)):
		if re.search(r"ttc_" + re.escape(str(t)) + r".*(\d+\.?\d*)$",file) == None:
			continue
		print file
		temperature=re.findall("ttc_(\d+\.?\d*)",file)[0]
		temps=re.findall("t(-?\d+\.?\d*)",file)[0]

		s=np.zeros((M,N))
		data=np.loadtxt(directory+file)
		for k in range(len(data)):     
			s[int(data[k,1]),int(data[k,0])]=int(data[k,2])


		plt.figure()
		plt.pcolormesh(s,rasterized=True,cmap="Pastel1") 
		plt.axis("image")
		plt.title("T="+temps)
		plt.savefig(directory + file + ".png",bbox_inches='tight')
		plt.close()

	os.system("convert -delay 50 -loop  $(ls " + directory + "*ttc_" + str(t) + "*.png | sort -tt -k4n,4) " + directory + "goutte_t" + str(t) + ".gif")

