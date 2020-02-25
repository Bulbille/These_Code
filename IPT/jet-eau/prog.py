from scipy import *
from scipy.optimize import*
from pylab import *

V1 = array([3,4,5,6,7,8,9,10,11])


V09 = array([3,4 ,5,6,7,8,9,10,11])
V08 = V09

V07 = array([3,4,5,6,7,8,9,10,11])
V06=V07

V05 = array([3,4,5,6,7,8,9,10,11])


A1  = array([100,100,51.069,51.069,40.262,35.362,23.811,13.397,9.888])
A09 = array([100,100,50.46,43.35,39.65,35.29,23.65,5.19,0.5])
A08 = array([100,100,49.720,42.793,31.482,24.048,8.473,3.026,0.5])
A07 = array([100,100,50.869,42.809,31.015,12.128,0.5,0.5,0.5])
A06 = array([100,100,35.480,19.283,10.560,3.977,0.5,0.5,0.5])
A05 = array([52.511,48.327,34.289,17.000,6.795,0.5,0.5,0.5,0.5])


R1  = exp(sqrt(A1/pi))
R09 = exp(sqrt(A09/pi))
R08 = exp(sqrt(A08/pi))
R07 = exp(sqrt(A07/pi))
R06 = exp(sqrt(A06/pi))
R05 = exp(sqrt(A05/pi))


Sep = np.array([1,0.9,0.8,0.7,0.6,0.5])


test = np.squeeze([A1,A09,A08,A07,A06,A05])
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


x = V1;y=Sep;

X, Y = np.meshgrid(x, y)
Z = test;


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z)

plt.show()






















plot(V1,R1,'o', color='blue',label="1 cm (exp)")
plot(V09,R09,'o', color='purple',label="0.9 cm (exp)")
plot(V09,R08,'o', color='black',label="0.8 cm (exp)")
plot(V07,R07,'o', color='grey',label="0.7 cm (exp)")
plot(V07,R06,'o', color='green',label="0.6 cm (exp)")
plot(V05,R05,'o', color='yellow',label="0.5 cm (exp)")


def Lineaire(x, a, b):
    return a*x+b

popt, pcov = curve_fit(Lineaire,V1 ,R1)
[a, b] = popt
plot(V1,a*V1+b,'-',color='blue' ,label="fit 1")

popt, pcov = curve_fit(Lineaire,V09 ,R09)
[a, b] = popt
plot(V09,a*V09+b,'-',color='purple' ,label="fit 0.9")

popt, pcov = curve_fit(Lineaire,V09 ,R08)
[a, b] = popt
plot(V09,a*V09+b,'-',color='black' ,label="fit 0.8")

popt, pcov = curve_fit(Lineaire,V07 ,R07)
[a, b] = popt
plot(V07,a*V07+b,'-',color='grey' ,label="fit 0.7")

popt, pcov = curve_fit(Lineaire,V07 ,R06)
[a, b] = popt
plot(V07,a*V07+b,'-',color='green' ,label="fit 0.6")

popt, pcov = curve_fit(Lineaire,V05 ,R05)
[a, b] = popt
plot(V05,a*V05+b,'-',color='yellow' ,label="fit 0.5")


title(" d en v ")
xlabel("V [kV]")
ylabel("exp(R) [cm]")
legend()
show()