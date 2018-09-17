#!/usr/bin/env python
# coding: utf-8

from pylab import *
import numpy as np
import matplotlib.pyplot as plt


R_Earth = 6378.137E3 # mètres
R_Jup   = 7.15E7 	 # mètres


f = 20 #Police des légendes des graphiques

data1 = genfromtxt("5_6.txt",dtype= float,skip_header=4,delimiter=None)
data2 = genfromtxt("10.txt",dtype= float,skip_header=4,delimiter=None)
data3 = genfromtxt("15.txt",dtype= float,skip_header=4,delimiter=None)

R1 = data1[:,0]*(R_Jup/R_Earth) #RAYON DE LA PLANÈTE CIBLE EN RAYON TERRESTRE
T1 = data1[:,1] #TEMPÉRATURE EFFECTIVE EN KELVIN
S1 = data1[:,2] #SÉPARATION ANGULAIRE EN BANDE 5.6µm 
F1 = data1[:,3]	# FLUX DES PLANÈTES EN BANDE 5.6µm

R2 = data2[:,0]*(R_Jup/R_Earth) #RAYON DE LA PLANÈTE CIBLE EN RAYON TERRESTRE
T2 = data2[:,1] #TEMPÉRATURE EFFECTIVE EN KELVIN
S2 = data2[:,2] #SÉPARATION ANGULAIRE EN BANDE 10µm 
F2 = data2[:,3]	# FLUX DES PLANÈTES EN BANDE 10µm

R3 = data3[:,0]*(R_Jup/R_Earth) #RAYON DE LA PLANÈTE CIBLE EN RAYON TERRESTRE
T3 = data3[:,1] #TEMPÉRATURE EFFECTIVE EN KELVIN
S3 = data3[:,2] #SÉPARATION ANGULAIRE EN BANDE 15µm 
F3 = data3[:,3]	# FLUX DES PLANÈTES EN BANDE 15µm

#F4 = data4a[:,0]

####################################################
#CALCUL DES HISTOGRAMMES
####################################################
#Rayon planétaire en rayon Terrestre

plt.title("Histogram of exoplanets for \n 5.6,10,15 $\mu$m",fontsize=f)
#plt.xlabel('Radius\n(Earth radius)',fontsize=f)
plt.xlabel('Temperature\n(Kelvin)',fontsize=f)
plt.ylabel('Number of planets',fontsize=f)

plt.title("Histogram of exoplanets for \n 5.6 $\mu$m",fontsize=f)
plt.xlabel('Temperature\n(Kelvin)',fontsize=f)
plt.ylabel('Number of planets',fontsize=f)
plt.hist(T1, bins='auto', color='blue',alpha =1)
grid(True)
plt.show()

plt.title("Histogram of exoplanets for \n 10 $\mu$m",fontsize=f)
plt.xlabel('Temperature\n(Kelvin)',fontsize=f)
plt.ylabel('Number of planets',fontsize=f)
plt.hist(T2, bins='auto', color='red',alpha = 1) 
grid(True)
plt.show()

plt.title("Histogram of exoplanets for \n 15 $\mu$m",fontsize=f)
plt.hist(T3, bins='auto', color='green',alpha = 1)
plt.xlabel('Temperature\n(Kelvin)',fontsize=f)
plt.ylabel('Number of planets',fontsize=f) 
grid(True)
plt.show()


####################################################
#FLUX PLANÉTAIRES VS SEPARATION ANGULAIRE
####################################################

figure(figsize=(15, 10), dpi=100)

title('PLANETARY FLUX \n vs \n ANGULAR SEPARATION',fontsize = f)
xlabel("Angular separation (mas)",fontsize = f)
ylabel("Flux (mJy)",fontsize = f)
xscale('log')
yscale('log')
grid(True)


plot(S1,F1,'bo',markersize=3,label="5.6$\mu$m")
plot(S2,F2,'ro',markersize=3,label="10$\mu$m")
plot(S3,F3,'ko',markersize=3,label="15$\mu$m")
	
legend(prop={'size':12})
	
#show()

#######################################################