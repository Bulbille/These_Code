#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=10
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate
from scipy.optimize import curve_fit
import numpy as np
import numpy.linalg as linalg
from pylab import *
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re

directory=sys.argv[1]
regex_chain="T-([\d.]+)(?:-\d)?(?:-t[\d.]*)?(?:-h-([\d.]+))?(?:-g-([\d.]+))?"
TC  = 2/np.log(1.+np.power(2,0.5))
J   = 1
LX  = 100
LY  = 20*2+1
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab','palegreen','lightseagreen','darkcyan','royalblue','navy','plum','coral','orange','olive','g','c','dodgerblue','violet','fuchsia','crimson','peru','goldenrod','forestgreen','aquamarine','aqua','blueviolet','pink','tomato','linen','antiquewhite','yellow','greenyellow','lime','cyan','indigo','darkmagenta']
fmt=["-",":","--", "-.","-",":","--", "-.","-",":","--", "-."] #linestyle
nb_color = size(colors)


###############################
# Liste des températures et 
# des longueur magnétiques 
###############################
temp = []
longueur = []
champ  = []
cisai  = []

for file in sorted(os.listdir(directory)):
    a=False
    b=False
    regex = re.search(regex_chain,file)
    if regex == None : 
        continue
    res = re.findall(regex_chain,file)
    if float(res[0][0]) not in temp:
        temp = np.append(temp,float(res[0][0]))
    if bool(res[0][1]) :
        if float(res[0][1]) not in champ :
            champ   = np.append(champ,float(res[0][1]))
        elif 0 not in champ:
            champ = np.append(champ,0.)
    if bool(res[0][2]) :
        if float(res[0][2]) not in cisai :
            cisai = np.append(cisai,float(res[0][2]))
    else :
        a=True;
        for i in cisai :
            if i <= 0.001:
                a=False
                break
        b=False
        for i in champ :
            if i<= 0.001:
                b=True
                break

##################################  :no
### Déclaration fonctions ########
##################################
zeros=sp.ai_zeros(3)[1]
def p(h,chi,zeros):
    nominateur = np.power(sp.airy(abs(h)/chi+zeros[0])[0],2)
    denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(abs(x)/chi+zeros[0])[0],2),0,np.inf)[0]
    return nominateur/denominateur


def calcul_moment(x,y):
    moment=[0,0,0]
    x_bin=abs(x[1]-x[0])
    x_bin=1.0
    moment[0] = np.sum(y)*x_bin
    for index, valeur in np.ndenumerate(y):
        index=index[0]
        moment[1]+= x[index]*valeur*x_bin
    moment[1] /= moment[0]
    for index, valeur in np.ndenumerate(y):
        index=index[0]
        moment[2]+=np.power(x[index]-moment[1],2)*valeur
    moment[2] = np.power(moment[2]*x_bin/moment[0],0.5)
    return moment


def translation(nb,m):
    return 1.*nb*m/(m-1)-1.*m/2

def transfert(kbt,champ,inter,l):
    d=np.zeros((l,l))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*champ*abs(translation(y1,l))/2 
                    -kbt*champ*abs(translation(y2,l))/2 
                    -kbt*inter*abs(translation(y1,l)-translation(y2,l)))
    return d

def fit_gauss(x,a,chi):
    return a*np.exp(-x**2/chi)

############################
# Moyenne des données
############################
big_data = {}

for file in sorted(os.listdir(directory)):
    if re.search(r"histo",file) != None :
        res = re.findall(regex_chain,file)	

        t=float(res[0][0]); f=cisai[0]; h=0
        if bool(res[0][1]):
            h=float(res[0][1])
        if bool(res[0][2]):
            f=float(res[0][2])
        try:
            data=np.loadtxt(directory+file)
            big_data[t,h,f] = data
        except: 
            continue
    else : 
        continue

def fit_exp(x,n,a):
    return a*np.exp(x/n)

###############################
# Accord TM et SOS
###############################

for j,t in enumerate(sorted(temp)) :
    BETA = 1/(t*TC)
    for k,h in enumerate(sorted(champ)):
        if h != 0.11:
            continue

        ########## Données Chipping #########
        courbe = [0,0]
        for g,f in enumerate(sorted(cisai)) :
            if f >= 2.56:
                continue
            y=[];num=[];
            try : 
                for point in big_data[t,h,f]:
                    if point[2] > 0 :
                       y=np.append(y,(point[0]+point[1])/2); num=np.append(num,point[2])
            except:
                continue
            moments = calcul_moment(y,num)
            print "Moments",t,h,f,moments
            y = (y-moments[1])/moments[2]; num  = np.multiply(num, moments[2]/moments[0])
            plt.plot(y,num,fmt[0], label="T ="+str(t)+", H ="+str(h)+", f ="+str(f),color=colors[g%nb_color])

#    lon = 2000
#    w, v    = linalg.eigh(transfert(BETA,0.33,J,2*int(lon)+1))
#    ytm   = np.arange(-int(lon),int(lon)+1,1)
#    tm  = np.power(v[:,np.argmax(w)],2)
#    moments = calcul_moment(ytm, tm)
#    ytm = (ytm-moments[1])/moments[2]; tm  = np.multiply(tm, moments[2]/moments[0])
#    plt.plot(ytm,tm,'-.',label='TM')
#    print moments

#y = np.linspace(-100,100,100000)
#airy = p(y,3,zeros)
#moments = calcul_moment(y, airy)
#y = (y-moments[1])/moments[2]; #airy  = np.multiply(airy, moments[2]/moments[0])
#plt.plot(y, airy,'black',label="Airy ",linewidth=2)

plt.ylabel('Height distribution')
plt.xlabel('Interface height)')
plt.yscale('log')
plt.xlim([-6,6])
plt.ylim([1e-7,5e-1])
plt.legend()
plt.savefig(directory+'kaw-ising2.pdf')
plt.show()
plt.close()


