#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fontsize=10
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import scipy as sy
from scipy.optimize import curve_fit
import numpy as np
from numpy import linalg as LA
import time
import sys
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab','palegreen','lightseagreen','darkcyan','royalblue','navy','plum','coral','orange','olive','g','c','dodgerblue','violet','fuchsia','crimson','peru','goldenrod','forestgreen','aquamarine','aqua','blueviolet','pink','tomato','linen','antiquewhite','yellow','greenyellow','lime','cyan','indigo','darkmagenta']


################################## 
### Déclaration fonctions ########
##################################

def powfit(x,a,n):
    return a*x**(-n)


############################
# Déclaration matrices

data=np.loadtxt(sys.argv[1])
#print data
print (data[-1,1])*pow(1/data[-1,0],-2)

plt.plot(data[:,0],data[:,1],'+',label="SOS F $= -\\frac{1}{\\beta}\ln(\lambda)-F(\infty)$",color="red")
plt.plot(data[:,0],(data[-1,1])*pow(data[:,0]/data[-1,0],-2),label="Power Law -2",color="red")


plt.xscale('log')
plt.yscale('log')
plt.axis([30,400,1e-4,1e-1])
plt.ylabel('$F(T,h,L)$')
plt.xlabel('Size $L_Y$')
plt.legend()
plt.show()
