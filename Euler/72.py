#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import numpy as np
from numpy import linalg as LA
import math 

N = 8

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

primes = np.concatenate(([1],primes),axis=0)

nbset = N-1
for n in np.arange(N+1) :
    for i,p in enumerate(primes) :
        if p >= n :
            break
        if n%p != 0 :
            nbset += n/p
            print n,p,i,n/p
print nbset
