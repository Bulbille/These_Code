#!/user/bin/python
# -*- coding: utf8

import matplotlib.pyplot as plt
import scipy as sy
import numpy as np
#from pylab import *
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import fnmatch
import re

y = np.arange(1,19)
ex = np.zeros(18)
ey = np.zeros(18)

nb_fichier=0
for file in sorted(os.listdir('./kaw75')): 
	if re.search(".+energie",file) == None: 
		continue
	nb_fichier+=1
	data=np.loadtxt('./kaw75/'+file)
#	print data[:,0]
	ex+=data[:,1]
	ey+=data[:,2]

ey/=nb_fichier
ex/=nb_fichier

plt.figure()
plt.subplot(211)
plt.plot(y,ex)

plt.subplot(221)
plt.plot(y,ey)

plt.show()
	
