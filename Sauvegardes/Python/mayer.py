import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from pylab import *
import scipy as sy
from scipy.interpolate import interp1d
from scipy.misc import derivative

### LECTURE DU FICHIER DES GRANDEURS THERMO###
## Lecture Glauber ##
with open('./glauber/thermo_kaw') as f:
        noms = f.readline().split('\t')
        dtipus = [('beta', sy.float32)] + [('mag', sy.float32)]+ [('energie', sy.float32)]+ [('chi', sy.float32)]+ [('c_v', sy.float32)]
        data = sy.loadtxt(f,delimiter=' ',dtype=dtipus)
g_beta=data['beta']
g_temp=1/g_beta
g_mag=data['mag']
g_energie=data['energie']
g_chi=data['chi']
g_c_v=data['c_v']
## Lecture kawasaki ##
with open('./kawasaki/thermo_kaw') as f:
        noms = f.readline().split('\t')
        dtipus = [('beta', sy.float32)] + [('mag', sy.float32)]+ [('energie', sy.float32)]+ [('chi', sy.float32)]+ [('c_v', sy.float32)]
        data = sy.loadtxt(f,delimiter=' ',dtype=dtipus)
k_beta = data['beta']
k_temp = 1/k_beta
k_mag = data['mag']
k_energie = data['energie']
k_chi = data['chi']
k_c_v = data['c_v']

## Interpolations
temp=np.linspace(max(g_temp[0],k_temp[0]), min(g_temp[-1],k_temp[-1]),80)

g_c_v = interp1d(g_temp,g_c_v, kind='linear')
g_chi = interp1d(g_temp,g_chi, kind='linear')
g_mag = interp1d(g_temp,g_mag, kind='linear', fill_value='extrapolate')

dmag = derivative(lambda temp: g_mag(temp),temp,dx=1e-6)


k_c_v = interp1d(k_temp,k_c_v, kind='linear')
k_chi = interp1d(k_temp,k_chi, kind='linear')
k_mag = interp1d(k_temp,k_mag, kind='linear', fill_value='extrapolate')

plt.figure(1)
plt.subplot(211)
plt.plot(temp,g_mag(temp),'o')
plt.plot(temp,dmag,'+')

plt.subplot(212)
gauche = g_c_v(temp)-k_c_v(temp)
droite = dmag*dmag*temp/g_chi(temp)
plt.plot(temp,gauche,'o')
plt.plot(temp,droite,'+')

plt.show()



#g_k_c = g_c_v[1:-1] - k_c_v[1:-1]
#partie_droite = g_temp[1:-1] / g_chi[1:-1] * dm_dt^2
