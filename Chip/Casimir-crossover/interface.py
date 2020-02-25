#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate
import numpy as np
import numpy.linalg as linalg
from pylab import *
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re

directory=sys.argv[1]
regex_chain="T-([\d.]+)(?:-\d)?(?:-t[\d.]*)?(?:-h-([\d.]+))?"
TC  = 2/np.log(1.+np.power(2,0.5))
J   = 1
LX  = 100
LY  = 20*2+1
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab','palegreen','lightseagreen','darkcyan','royalblue','navy','plum','coral','orange','olive','g','c','dodgerblue','violet','fuchsia','crimson','peru','goldenrod','forestgreen','aquamarine','aqua','blueviolet','pink','tomato','lighsalmon','linen','antiquewhite','yellow','greenyellow','lime','cyan','indigo','darkmagenta']
fmt=["-",":","--", "-.","-",":","--", "-.","-",":","--", "-."] #linestyle


###############################
# Liste des températures et 
# des champs magnétiques 
###############################
temp = []
champs = []

for file in sorted(os.listdir(directory)):
	a=False
	regex = re.search(regex_chain,file)
	if regex == None : 
		continue
	res = re.findall(regex_chain,file)
	if float(res[0][0]) not in temp:
		temp = np.append(temp,float(res[0][0]))
	if  bool(res[0][1]) :
		if float(res[0][1]) not in champs :
			champs = np.append(champs,float(res[0][1]))
	else :
		a=True;
		for i in champs :
			if i<= 0.001:
				a=False
				break
	if a :
		champs = np.append(champs,0.)

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


############################
# Moyenne des données
############################
big_data = {}
nb = {}

for file in sorted(os.listdir(directory)):
	if re.search(r"histo",file) != None :
		res = re.findall(regex_chain,file)	

		t=float(res[0][0]); h=0
		if bool(res[0][1]):
			h=float(res[0][1])
		data=np.loadtxt(directory+file)
		data=data[:,:][~np.isnan(data[:,1])]
		if (t,h) not in big_data:
			big_data[t,h] = data
			nb[t,h]=1
		else :
			big_data[t,h] = np.add(data,big_data[t,h])
			nb[t,h] += 1
	else : 
		continue

for points in big_data :
	big_data[points] /= nb[points]

###############################
# Calcul des moments
###############################
zeros=sp.ai_zeros(3)
chi=np.power(integrate.quad(lambda x : np.power(sp.airy(abs(x)+zeros[1][0])[0],2),0,np.inf)[0]/integrate.quad(lambda x : np.power(x*sp.airy(abs(x)+zeros[1][0])[0],2),0,np.inf)[0],0.5)

for t in temp :
    BETA = 1/(t*TC)
    for i,h in enumerate(champs) :
        if h != 0.1 :
            continue
        y=[];valeur=[];
        for point in big_data[t,h]:
            if abs(point[1] > 1):
                y=np.append(y,point[0]); valeur=np.append(valeur,point[1])
        moments = calcul_moment(y,valeur)
        valeur = np.multiply(valeur, 1/moments[0])
       # x=(x-moments[1])/moments[2]; valeur = np.multiply(valeur, moments[2]/moments[0])
       # valeur = np.multiply(valeur, 1/moments[0])
        plt.plot(y,valeur,label="t="+str(t)+",h="+str(h),color=colors[i])
        print h, moments[2]

        w, v    = linalg.eigh(transfert(BETA,h,J,LY))
        ################# rajouter les sommes avec LX !!
        valeur  = sum( np.power(v[:,index],2)*np.power(l,LX) for index,l in enumerate(w))
        norma   = sum( np.power(l,LX) for l in w)
        valeur  = np.multiply(valeur,1/norma)
        moments = calcul_moment(range(LY),valeur)
        y       = np.arange(LY);             moments = calcul_moment(y,valeur)
        y       = y-moments[1]
       # x       = (x-moments[1])/moments[2]; valeur  = np.multiply(valeur, moments[2]/moments[0])
        plt.plot(y,valeur,'-.', label="t="+str(t)+",h="+str(h),color=colors[i])
        break

#xairy=np.linspace(min(x),max(x),50)
#plt.plot(xairy,p(xairy,chi,zeros[1]),".", label="airy", color="black")
#plt.plot(xairy,np.exp(-xairy**2/2)/np.power(2*np.pi,0.5),"-+",label="gaussienne",color="black")
plt.axis([min(y),max(y),min(valeur),max(valeur)])
plt.yscale('log')
plt.legend()
plt.show()
plt.close()
