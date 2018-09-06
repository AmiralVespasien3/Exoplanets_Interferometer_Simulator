#!/usr/bin/env python
# coding: utf-8
# programme.py
# variables.py

from matplotlib.pyplot import *
from variables import *
from fonctions import C_sun

#############################################################

def Rayon_Temperature(f,C,R):
	figure(figsize=(15, 10), dpi=100)
	title(u"Radius vs Temperature ",fontsize = f)
	xlabel(u"Effective temperature (°K)",fontsize = f)
	ylabel(u"Planetary radius (RJ)",fontsize = f)
	xscale('log')
	yscale('log')
	#xlim(xmin=0)
	#ylim(ymin=1E-17)
	grid(True)

	plot(T_pr,R_pr/RJ,'bo',markersize=2)#A modifier

	T1 = []; R1 = []; T2 = []; R2 = []

	for i in range(0,len(C)):
		#print C[i]
		if 1e-5<=C[i]<= 1e-4:
			T1.append(T_pr[i])
			R1.append(R_pr[i]/RJ)
		if C[i]> 1e-5:
			T2.append(T_pr[i])
			R2.append(R_pr[i]/RJ)

	plot(T1,R1,'ro',markersize=4,label="1e-5<Contrast<1e-4",alpha=1)
	plot(T2,R2, 'go',markersize=3,label="Contrast>1e-5",alpha=1)

	#M0V
	R_rouge=6.2*RJ
	T_rouge=4e3


	#M9V
	R_rouge=0.8*RJ
	T_rouge=2e3


	T_lim = pow(  ( pow(R_rouge,2)*(1-A_B) ) / (4*(R*m2r)**2)  ,   1./4  ) * T_rouge 
	#axvline(x=T_lim,color='red')

	S_plim = (R_rouge/2)*(1-A_B)**(1./2) * (T_rouge/T_lim)**2
	#print "S_plim  " + str(S_plim/m2r)

	legend(prop={'size':12})
	show()

###############################################################

def Hemisphere(C,f,**data):

	figure(figsize=(15, 10), dpi=100)

	#TEST HEMISPHERE

	title('Contrast vs Angular separation \n'+ "$\lambda$ = " + str("%.2e"%L) + "m",fontsize = f)
	xlabel("Angular separation (mas)",fontsize = f)
	ylabel("Contrast (None)",fontsize = f)
	xscale('log')
	yscale('log')
	ylim(ymin=1E-17)
	grid(True)

	for i in range(0,len(P_d)):
		if P_d[i] > 0:
			data["S_pa_N"].append(S_pa[i])
			data["C_N"].append(C[i])
			
		else :
			data["S_pa_S"].append(S_pa[i])
			data["C_S"].append(C[i])

	plot(data["S_pa_N"],data["C_N"],'bo',markersize=3,label="North Hemisphere visible planets")
	plot(data["S_pa_S"],data["C_S"],'ro',markersize=3,label="South Hemisphere visible planets")

	plot(100,C_sun(L,R_Terre,T_Terre),'gs',markersize=5,label="Earth at 10pc")
	plot(500,C_sun(L,RJ,T_Terre),'k*',markersize=5,label="Jupiter at 10pc")
	legend(prop={'size':12})
	show()
	

################################################################

def Proxima(PROXIMA_R,PROXIMA_T,f):


	R_Terre = 6378.137E3 # mètres
	M_Terre = 5.97E24 # kg
	R_lim = RJ #Ryon de Jupiter
	T_lim = T_pr[3419] # WASP-18b  2000°K : i = 3419

	print RJ/R_Terre
	print T_lim

	#Initialisation des listes
	PROXIMA_R.append(R_pr[3338])
	PROXIMA_T.append(T_pr[3338])

	i = 0
	r = 0.1*R_Terre
	t = 10

	while PROXIMA_R[i] <= RJ and PROXIMA_T[i] <= T_lim :
		PROXIMA_R.append(PROXIMA_R[i]+r)
		PROXIMA_T.append(PROXIMA_T[i]+t)  
		i += 1

	# Convertion des rayons en rayons Terrestres
	for i in range(0,len(PROXIMA_R)):
		PROXIMA_R[i] = PROXIMA_R[i]/R_Terre
		#print str("%.4e"%PROXIMA_R[i]) + '\t' + str("%.4e"%PROXIMA_T[i]) 
		#print str("%3.e"%PROXIMA_R[i])

	#GRAPHIQUE
	figure(figsize=(15, 10), dpi=100)
	title(u"Radius vs Temperature ",fontsize = f)
	xlabel(u"Effective temperature (°K)",fontsize = f)
	ylabel(u"Planetary radius (R_Earth)",fontsize = f)
	#xscale('log')
	#yscale('log')
	#xlim(xmin=PROXIMA_T[-1])
	#ylim(ymin=12)
	grid(True)
	
###############################################################

def Paranal(C,f,R,**data):

	figure(figsize=(15, 10), dpi=100)

	##Constraste vs Separation Angulaire
	title("Contrast vs Angular separation\n" + "$\lambda$ = " + str("%.2e"%L) + " m",fontsize = f)
	xlabel("Angular separation (mas)",fontsize = f)
	ylabel("Contrast (None)",fontsize = f)
	xscale('log')
	yscale('log')
	ylim(ymin=1E-17)
	grid(True)

	#TEST PARANAL

	for i in range(0,len(P_d)):
		if P_d[i] <= 0:
			data["S_pa_S"].append(S_pa[i])
			data["C_S"].append(C[i])

	for i in range(0,len(P_d)):
		if -84 <= P_d[i] <= 36 :
			data["S_pa_P"].append(S_pa[i])
			data["C_P"].append(C[i])

	
	plot(data["S_pa_P"],data["C_P"],'bo',markersize=3,label="Visible planets from Paranal")
	axvline(x=R,color='black')
	axhline(y=C_lim,color='black')
	axhline(y=C_lim/10,color='black')

	plot(100,1E-10,'gs',markersize=5,label="Earth at 10pc")
	plot(500,3.9E-17,'k*',markersize=5,label="Jupiter at 10pc")

	legend(prop={'size':12})
	show()

###############################################################

def Flux_Distance(D_s,F_p):
	f = 20
	figure(figsize=(8, 6), dpi=80)

	title(u"Flux des planètes vs éloignement système",fontsize = f)
	xlabel(u"Eloignement du système (pc)",fontsize = f)
	ylabel("Flux de la pLanete (Jy)",fontsize = f)
	xscale('log')
	yscale('log')
	ylim(ymin=1E-17)
	grid(True)
	plot(D_s/pc,F_p,'ko',markersize=3)
	plot(1,1.19E-9,'bo',markersize=4)
	plot(5,4.57E-16,'ro',markersize=4)

	show()

###############################################################