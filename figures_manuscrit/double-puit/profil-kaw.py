#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fontsize=14
plt.rc("font",size=fontsize)

puit = plt.subplot(111)
import numpy as np

def phi4(tab,i,dh):
    n = i-2
    a = 12*tab[n]**2*(tab[n+1]-tab[n-1])
    b = 4*tab[n-1]-6*tab[n]+4*tab[n+1]-tab[n-2]
    return float(a+b-4*dh)


L = 100
dx = 0.0001
fun = np.empty(L)
fun[0] = fun[1] = fun[2 ] = fun[3] = 0

for x in np.arange(L)[4:]:
    fun[x] = phi4(fun,x,dx)
    print(x,fun[x])

puit.plot(np.arange(L),fun)

