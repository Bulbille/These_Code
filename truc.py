
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def PFD(V,t):
    vx,vy = V
dvdx = -B*vy/np.sqrt(vx**2+vy**2 +e) - C*vx
    dvdy = -A + B*vx/np.sqrt(vx**2+vy**2+e) - C*vy
        return(dvdx,dvdy)

TMIN = 0
TMAX = 2
TN   = 20000
e    = 1e-5
A,B,C= (1,0,1)
V    = [1,0.5]
t    = np.linspace(TMIN,TMAX,TN)
v    = odeint(PFD,V,t)
vx   = v[:,0]
vy   = v[:,1]

def integration(v,tmin,tmax,tn):
    x = np.empty(int(tn/2))
for i,zz in enumerate(x) :
    tmp = 0
          plt.show()

