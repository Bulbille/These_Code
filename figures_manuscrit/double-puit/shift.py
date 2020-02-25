#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fontsize=14
plt.rc("font",size=fontsize)

puit = plt.subplot(111)
import numpy as np

def doublepuit(x):
    return (x**2-1)**2
def lin(x,a):
    return a*x

L = 1.5
x = np.linspace(-L,L,400)
puit.plot(x,doublepuit(x))
puit.plot(x,doublepuit(x)+lin(x,1))

puit.set_xlabel('$\phi$')
puit.set_ylabel('$V(\phi)$')

plt.savefig('shift.pdf')
plt.show()
