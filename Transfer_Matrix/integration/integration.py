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

def sos(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*champ*(abs(L-y1)+abs(L-y2))/2
                    -kbt*inter*abs(y1-y2) )
    return d

def gauss(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*champ*(abs(L-y1)+abs(L-y2))/2
                    -kbt*inter*pow(y1-y2,2) )
    return d

def matrice_phi(L):
    d=np.zeros((2*L+1,2*L+1))
    for y,i in enumerate(d):
        d[y][y] = abs(L-y)
    return d

def powfit(x,a,n,b):
    return a*x**(-n)+b


############################
# Déclaration matrices

#data=np.loadtxt(sys.argv[1])
#plt.plot(data[:,1],data[:,2],label="Glauber")
#data=np.loadtxt(sys.argv[2])
#plt.plot(data[:,1],data[:,2],label="Kawasaki")

TC  = 2/np.log(1.+np.power(2,0.5)); BETA = 1/TC ; J = 1; 
R = np.exp(-BETA); epsilon = (1+R)/(1-R); maxsos = -np.log(epsilon)/BETA
epsgauss = 1+2*np.sum(pow(R,pow(np.arange(1,10),2))) ; maxgauss = -np.log(epsgauss)/BETA

Hmin = 0.0; Hn = 2000; Hmax = 30; 
champs = np.linspace(Hmin,Hmax,Hn)

LX = 50
LY = 10
lys = np.arange(30,60)

matrix  = np.empty([np.size(lys),2])
esos = np.empty([np.size(lys),2])
egauss = np.empty([np.size(lys),2])

#for k,h in enumerate(champs):
for l,LY in enumerate(lys) :
    print LY
    tm = sos(BETA,Hmin,J,LY)
    w,v = LA.eigh(tm)

    esos[l] = [LY,-1/BETA*np.log(max(w))]

    tm = gauss(BETA,0,J,LY)
    w,v = LA.eigh(tm)

    egauss[l] = [LY,-1/BETA*np.log(max(w))]

#    omega = 0
#    matphi = matrice_phi(LY)
#    for k,h in enumerate(np.linspace(Hmin,Hmax,Hn)):
#        tm = sos(BETA,h,J,LY)
#        w,v = LA.eigh(tm)
#        Z = sum(w**LX)
#        tmpow = LA.matrix_power(tm,LX)
#        trace = np.trace(np.dot(matphi,tmpow))
#        phi = trace/Z
#        if(k == 0 or k == Hn-1):
#            omega += phi
#        elif (k%2 == 0 ):
#            omega += 2*phi
#        else:
#            omega += 4*phi
#    omega /= (3*Hn)/(Hmax-Hmin)
#    matrix[l] = [LY,-omega]
 #   Z = sum(w**LX)
 #   tmpow = LA.matrix_power(tm,LX)
 #   trace = np.trace(np.dot(matphi,tmpow))

 #   matrix[l] = [LY,trace/Z]

#plt.plot(matrix[:,0],matrix[:,1],label="I $=-\\int_0^\infty M(B) dB$")
#popt,pcov = curve_fit(powfit,matrix[:,0],matrix[:,1],p0=[1,2,-4])
#plt.plot(matrix[:,0],powfit(matrix[:,0],popt[0],popt[1],popt[2]),label="Fit I $="+str(round(popt[0],3))+"\\times LY^{-"+str(round(popt[1],3))+"}"+str(round(popt[2],3))+"$")
#

plt.plot(esos[:,0],esos[:,1]-maxsos,'+',label="SOS F $= -\\frac{1}{\\beta}\ln(\lambda)-F(\infty)$",color="red")
plt.plot(esos[:,0],(esos[-1,1]-maxsos)*pow(esos[:,0]/esos[-1,0],-2),label="Power Law -2",color="red")
#popt,pcov = curve_fit(powfit,esos[:,0],esos[:,1],p0=[1,2,-4])
#plt.plot(esos[:,0],powfit(esos[:,0],popt[0],popt[1],popt[2]),'+',label="Fit F $="+str(round(popt[0],3))+"\\times LY^{-"+str(round(popt[1],3))+"}"+str(round(popt[2],3))+"$")

plt.plot(egauss[:,0],egauss[:,1]-maxgauss,'+',label="Gauss F $= -\\frac{1}{\\beta}\ln(\lambda)$",color="Blue")
plt.plot(egauss[:,0],(egauss[-1,1]-maxgauss)*pow(egauss[:,0]/egauss[-1,0],-2),label="Power Law -2",color="Blue")
#popt,pcov = curve_fit(powfit,egauss[:,0],egauss[:,1],p0=[1,2,-4])
#plt.plot(egauss[:,0],powfit(egauss[:,0],popt[0],popt[1],popt[2]),'+',label="Fit F $="+str(round(popt[0],3))+"\\times LY^{-"+str(round(popt[1],3))+"}"+str(round(popt[2],3))+"$")

np.savetxt('sos',esos)
np.savetxt('gauss',egauss)
esos[:,1] = esos[:,1]-maxsos
egauss[:,1] = egauss[:,1]-maxgauss
np.savetxt('sos-offset',esos)
np.savetxt('gauss-offset',egauss)

#print matrix,energie

plt.ylabel('$F(T,h,L)$')
plt.xlabel('Size $L_Y$')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.savefig('gaussian-glauber.pdf')
plt.show()
plt.close()
