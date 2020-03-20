#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fontsize=8
plt.rc("font",size=fontsize)

puit = plt.subplot(111)
import numpy as np

def doublepuit(x,m,l):
    return m*x**2+l**4/24*x**4

L = 6
x = np.linspace(-L,L,400)

for m in np.round(np.linspace(-1,0.1,5),2) :
    print(m)
    puit.plot(x,doublepuit(x,m,1),label="$T-T_C="+str(m)+"$")
puit.plot(x,doublepuit(x,-1,1)+1.5*x,'-.',color="black",label="$T-T_C=-1.0$, avec champ magnétique")

puit.legend()
puit.set_xlabel('$\phi$')
puit.set_ylabel('$V(\phi)$')

plt.savefig('double-puit-en-fonction-temp.pdf')
plt.show()
