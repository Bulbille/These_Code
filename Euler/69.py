#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import numpy as np
from numpy import linalg as LA
import math 

N = 1e6

## D'abord il faut trouver tous les nombres premiers
numbers = np.arange(N+1)
p = 2
while p < N**0.5 :
    if numbers[p] != 0 :
        tmp = 2
        while p*tmp <= N :
            numbers[p*tmp] = 0
            tmp +=1
    p += 1
primes = np.where(numbers != 0)[0][1:]
#### Utilisationd des propriétés connues pour préremplir le tableau
indicatrice = np.arange(N+1)

for p in primes :
    tmp = 1
    while p*tmp <= N :
        indicatrice[p*tmp] *= 1.-1./p
        tmp += 1


for i,j in enumerate(indicatrice) :
    if i == 0 :
        pass
    indicatrice[i] = i/j
print np.argmax(indicatrice[2:])+2
exit()

import time
temps = time.clock()
## phi(p^k) = p^k - p^k-1
for p in primes :
    k = 1
    while pow(p,k) <= N :
        indicatrice[pow(p,k)] = pow(p,k)-pow(p,k-1)
        k += 1
## phi(p*p') = phi(p)*phi(p')
for i,p in enumerate(primes) : 
    for p2 in primes[i+1:] :
        if p*p2 > N :
            break
        indicatrice[p*p2] = indicatrice[p]*indicatrice[p2]

for i,p in enumerate(primes) :
    k = 1
    isOkay = True
    while isOkay :
#        print k
        for p2 in primes[i+1:] :
            nb = pow(p,k)*p2
            if nb > N :
                isOkay = False
                break
            indicatrice[nb] = nb*(1-p)*(1-p2)
        k += 1

print indicatrice
exit()

## Fonction totient
def totient(x) :
    phi = x
    if x in primes :
        return x/(x-1)
    for p in primes : 
        if p > x/2 :
            break
        if x%p ==0 :
            phi *= (1-1./p)
    if x in primes : 
        phi *= (1-1/x)
    return phi
## Calcul pour tout le reste
for i,j in enumerate(indicatrice) :
    if j == 0 :
        indicatrice[i] = totient(i)
for i,j in enumerate(indicatrice) :
    indicatrice[i] = i/j

print np.argmax(indicatrice[2:])+2
