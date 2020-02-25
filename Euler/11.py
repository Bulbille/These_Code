#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import numpy as np
from numpy import linalg as LA

grid = np.loadtxt('nbproj11')
taille = len(grid)
print grid
print grid[0][0],grid[taille-1][taille-1]

LARG = 4
maxi = 0

#Vertical
for x in np.arange(taille) :
    for y in np.arange(taille-LARG) :
        tmp = grid[x][y]
        for i in np.arange(1,LARG):
            tmp *= grid[x][y+i]
        if tmp > maxi :
            maxi = tmp
print maxi
#horitonztal
for y in np.arange(taille) :
    for x in np.arange(taille-LARG) :
        tmp = grid[x][y]
        for i in np.arange(1,LARG):
            tmp *= grid[x+i][y]
        if tmp > maxi :
            maxi = tmp
print maxi

#diagonales

#gauche
for y in np.arange(taille):
    for x in np.arange(taille-LARG) :
        if y>taille-LARG or x > taille-LARG :
            continue
        tmp = grid[x][y]
        for i in np.arange(1,LARG):
            tmp *= grid[x+i][y+i]
        if tmp > maxi :
            maxi = tmp
print maxi
#droite
for y in np.arange(taille):
    for x in np.arange(taille-LARG) :
        if y>taille-LARG or x > taille-LARG :
            continue
        tmp = grid[x][y]
        for i in np.arange(1,LARG):
            tmp *= grid[x-i][y+i]
        if tmp > maxi :
            maxi = tmp
print int(maxi)
