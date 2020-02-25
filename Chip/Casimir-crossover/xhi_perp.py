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
regex_chain="T-([\d.]+)(?:-\d)?(?:-t[\d.]*)?(?:-h-([\d.]+))?(?:-lx-([\d.]+))?"
TC  = 2/np.log(1.+np.power(2,0.5))
J   = 1
LX  = 100
LY  = 20*2+1
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab','palegreen','lightseagreen','darkcyan','royalblue','navy','plum','coral','orange','olive','g','c','dodgerblue','violet','fuchsia','crimson','peru','goldenrod','forestgreen','aquamarine','aqua','blueviolet','pink','tomato','lighsalmon','linen','antiquewhite','yellow','greenyellow','lime','cyan','indigo','darkmagenta']
fmt=["-",":","--", "-.","-",":","--", "-.","-",":","--", "-."] #linestyle


###############################
# Liste des températures et 
# des longueur magnétiques 
###############################
temp = []
longueur = []
champ  = []

for file in sorted(os.listdir(directory)):
    a=False
    b=False
    regex = re.search(regex_chain,file)
    if regex == None : 
        continue
    res = re.findall(regex_chain,file)
    if float(res[0][0]) not in temp:
        temp = np.append(temp,float(res[0][0]))
    if  bool(res[0][1]) :
        if float(res[0][1]) not in champ :
            champ   = np.append(champ,float(res[0][1]))
    if  bool(res[0][2]) :
        if float(res[0][2]) not in longueur :
            longueur = np.append(longueur,float(res[0][2]))
    else :
        a=True;
        for i in longueur :
            if i<= 0.001:
                a=False
                break
        b=True
        for i in champ :
            if i<= 0.001:
                b=False
                break
    if a :
        longueur = np.append(longueur,0.)
    if b :
        champ = np.append(champ,0.)

##################################  :no
### Déclaration fonctions ########
##################################

def p(h,chi,zeros):
        nominateur = np.power(sp.airy(abs(h)/chi+zeros[0])[0],2)
        denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(abs(x)/chi+zeros[0])[0],2),0,np.inf)[0]
        return nominateur/denominateur
def calcul_moment(x,y):
    moment=[]
    x_bin=abs(x[1]-x[0])
    moment=np.append(moment,np.sum(y)*x_bin)
    moment=np.append(moment,0)
    for index, valeur in np.ndenumerate(y):
        index=index[0]
        moment[1]+= x[index]*valeur*x_bin
    moment[1] /= moment[0]

    moment=np.append(moment,0)
    for index, valeur in np.ndenumerate(y):
        index=index[0]
        moment[2]+=np.power(x[index]-moment[1],2)*valeur*x_bin
    moment[2] = np.power(moment[2]/moment[0],0.5)

    moment=np.append(moment,0)
    for index, valeur in np.ndenumerate(y):
        index=index[0]
        moment[3]+=np.power((x[index]-moment[1])/moment[2],3)*valeur*x_bin
    moment[3] /= moment[0]
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

def fit_pow(x,a,n):
    return a*np.power(x,n)

def largeur_mi_hauteur(x,y):
    miH = max(y)/2
    #Xmoins
    xm   = min(np.where(num > max(num)/2)[0])
    yxm  = y[xm]
    yxmp = y[xm-1]
    a    = (yxm-yxmp)
    b    = yxm - a*xm
    xmoins  = (miH-b)/a
    #Xplus
    xp   = max(np.where(num > max(num)/2)[0])
    yxp  = y[xp]
    yxpp = y[xp-1]
    a    = (yxp-yxpp)
    b    = yxp - a*xp
    xplus  = (miH-b)/a

    return (xplus-xmoins)/2

############################
# Moyenne des données
############################
big_data = {}

for file in sorted(os.listdir(directory)):
    regex = re.search(regex_chain,file)
    if regex == None : 
        continue
    res = re.findall(regex_chain,file)
    temp    = float(res[0][0])
    mu      = float(res[0][1])
    if mu not in muspace:
        muspace = np.append(muspace,mu)
    if temp not in tempspace:
        tempspace = np.append(tempspace,temp)
    if re.search('histo',file) :
        big_histo[temp,mu] = np.loadtxt(directory+file)
    elif re.search('snap',file) :
        big_snap[temp,mu] = np.loadtxt(directory+file)


for file in sorted(os.listdir(directory)):
    if re.search(r"histo",file) != None :
        res = re.findall(regex_chain,file)	

    t=float(res[0][0]); l=0; h=0
    if bool(res[0][1]):
        h=float(res[0][1])
    if bool(res[0][2]):
        l=float(res[0][2])
        try:
            data=np.loadtxt(directory+file)
            data=data[:,:][~np.isnan(data[:,1])]
            big_data[t,h,l] = data
        except: 
            continue
    else : 
        continue

###############################
# Calcul des moments
###############################
zeros=sp.ai_zeros(3)
chi=np.power(integrate.quad(lambda x : np.power(sp.airy(abs(x)+zeros[1][0])[0],2),0,np.inf)[0]/integrate.quad(lambda x : np.power(x*sp.airy(abs(x)+zeros[1][0])[0],2),0,np.inf)[0],0.5)

for j,t in enumerate(temp) :
    yt      = np.arange(LY) - LY/2
    BETA = 1/(t*TC)
    for k,h in enumerate(champ):
        w, v    = linalg.eigh(transfert(BETA,h,J,LY))
        yt       = yt-LY/2

        courbe = [0,0]

        for i,lon in enumerate(longueur) :
            y=[];num=[];
            try : 
                for point in big_data[t,h,lon]:
                    y=np.append(y,point[0]); num=np.append(num,point[1])
            except:
                continue

            tm  = sum( np.power(v[:,index],2)*np.power(lam,lon) for index,lam in enumerate(w)) / sum( np.power(lam,lon) for lam in w)
            val = abs( largeur_mi_hauteur(y,tm) - largeur_mi_hauteur(y,num))
#            courbe = np.vstack([courbe,[lon,largeur_mi_hauteur(y,num)]])
            courbe = np.vstack([courbe,[lon,val]])

        courbe = np.delete(courbe,0,0)
        courbe = courbe[courbe[:,0].argsort()]

        plt.plot(courbe[:,0],courbe[:,1],color=colors[k],label="t="+str(t)+"h="+str(h))

        rows = np.where(courbe[:,0] > 10)
        popt,pcov = curve_fit(fit_pow,courbe[rows][:,0],courbe[rows][:,1],p0=[1,-1])
        plt.plot(courbe[:,0],fit_pow(courbe[:,0],popt[0],popt[1]),'-.',color=colors[k],label="Exponent : "+str(+popt[1]))

        print h,popt

#plt.axis([min(courbe[:,0]),max(courbe[:,0]),min(courbe[:,1]),max(courbe[:,1])])
plt.ylabel('Largeur mi-hauteur de p(h)')
plt.xlabel('Size system LX (LY = 20)')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig('size-effectsmiH.pdf')
plt.show()
plt.close()
