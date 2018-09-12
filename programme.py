#!/usr/bin/env python
# coding: utf-8

import warnings
import sys

import encodings
import numpy as np

import time
start_time = time.time()

from pylab import *
from variables import *
from fonctions import *
from Graphique import *


if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

###############################################################################
#CALCUL DES FLUX
###############################################################################
# THERMAL FLUX
F_s = F_s(L,*T_s)
F_p = F_pa(L,Nr_p,T_pr,R_pr,D_s)
# REFLECTIVE FLUX
F_rp  = A_B * pi * (R_s/D_s)**2 * F_s * (R_pr**2)/(4*O_p**2)

###############################################################################
#ECRITURE DES FLUX DANS UN FICHIER DE SORTIE
###############################################################################
fichier = open("./Intefero_spatial/planets_data.tbl", "w")
fichier.write('This file contains all the planetatry fluxes in may 2018 \n')
fichier.write('Flux(mJy)\n')
for i in range(0,len(F_p)):
	fichier.write(str("%.3e" %( F_p[i]/(1e-3)) ) +'\n')
fichier.close()


###############################################################################

###############################################################################
#CALCUL DU CONTRASTE 
###############################################################################

C = (F_p+F_rp)/F_s

##############################################################################
#CATALOGUE VIKING########

def title_label(a,b,c,f):
	figure(figsize=(15, 10), dpi=100)
	title(a,fontsize = f)
	xlabel(b,fontsize = f)
	ylabel(c,fontsize = f)
def log_label():
	#xscale('log')
	yscale('log')
	grid(True)


''' Download from Viking targets ! '''
a = ['HD179949', 'HD75289', 'HD102195', 'HD77338', 'HD285507', 'GJ86', 'HD168746', 'HD10180', 'GJ876', 'HD27894', '61Vir', 'GJ3998', 'HD1461', 'HD40307', 'HD160691', 'nuOph', 'HD168443', 'HD192263', 'HD102117', 'HD16417']
b = ['b',			'b'			,'b'		,'b'		,'b'	,'b'		,'b'        ,'c'	  ,'d'	     ,'b'     ,'b'    ,'b'      ,'b'        ,'b'      ,'c'     ,'b'     ,'b'         ,'b'         ,'b'         ,'b']

#####################################################################################################################################
'''
figure(figsize=(15, 10), dpi=100)

f =20

title(u"VIKING TARGETS \\ Constrast vs Angular separation ",fontsize = f)
xlabel(u"Angular separartion (mas)",fontsize = f)
ylabel(u"Contrast",fontsize = f)

xscale('log')
yscale('log')
grid(True)

CV = []; SV = []
for i in range(0,len(a)):
	p = planetes_names(a[i],b[i])
	#print '\n'
	#print a[i],'\t',b[i],'\t',p
	#print 'Mass','\t','Temp','\t\t','Cont','\t\t','Sep'
	#print M_p[p]/MJ,'\t',"%.3e" %T_pr[p],'\t',"%.3e" %C[p],'\t',"%.3e" %S_pa[p]
	CV.append(C[p])
	SV.append(S_pa[p])
	plot(S_pa[p],C[p],'o',label= a[i]+'\t'+b[i])

figure(figsize=(15, 10), dpi=100)

t1 = u"VIKING TARGETS \\ Mass vs Effective Temperature "
t2 = u"Effective Temperature (°K)"
t3 = u"Mass (Jovian mass)"

title_label(t1,t2,t3,f)
log_label()


for i in range(0,len(a)):
	p = planetes_names(a[i],b[i])
	plot(T_pr[p],M_p[p]/MJ,'o',label= a[i]+'\t'+b[i])

legend(prop={'size':10})
#show()
'''
###############################################################################
# TRIAGE DES ARRAYS SELON LA REFERENCE DE L'ARRAY DE CONTRASTE TRIÉ DANS L'ORDRE
# CROISSANT
###############################################################################

'''
sort2L(C,F_p)
print("--- %s seconds ---" % (time.time() - start_time))
'''


###############################################################################
# CALCUL DU TEMPS POUR INTÉGRER LA SOURCE (EN JOURNÉE)
###############################################################################
#S = 1e-15 #Sensibilité des Australiens

Time = np.array( zeros( len(Nr_p) ) )
for i in range(0,len(Time)):
	Time[i] = S/F_p[i]
	Time[i] /= Trans
#sort2L(C,Time)



################################################################################
#ALGORTITHME DE RECHERCHE DE CIBLES POUR Hi5
################################################################################
Cibles(L,R,C,F_p,Time,T_max)

#################################################################################
#GRAPHIQUE
#################################################################################

#Police des axes
f = 20 

figure(figsize=(15, 10), dpi=100)
title(u"CONTRAST MAP " + "\n$\lambda$ ="+ str("%.2e"%L)+"m" ,fontsize = f)
xlabel(u"Effective temperature (°K)",fontsize = f)
ylabel(u"Planetary radius (R_Earth)",fontsize = f)
xscale('log')
yscale('log')
grid(True)

C2 = []
T_pr2 = []
R_pr2 = []

for i in range(0,len(C)):
	if T_pr[i] > 90 :
		C2.append(C[i])
		T_pr2.append(T_pr[i])
		R_pr2.append(R_pr[i]/R_Terre)

array(C2)
array(T_pr2)
array(R_pr2)

scatter(T_pr2,R_pr2,c=log(C2),s=75,alpha=0.5)  # Traçage des points représentant les exoplanètes connues
legend(prop={'size':f})
cbar = colorbar()
cbar.set_label('log10(C)', rotation=270, fontsize = f)

#show()

###########################################################################
#Rayon de planètes simulées vs Températures effectives 
#Sur la base des don* (R_pr[i]/D_pr[i])**2 *pinées de Proxima Centauri b

#Rayon des plnètes vs Températures effectives planètaires
#PROXIMA CENT,AURI b : i = 3338
#Selon Fresin : R_super_Terre entre 1.4 et 2 R_terrestre

#############################################################################################
#CALCULS DES COURBES D'ISOCONTRASTES : C=k
##############################################################################################
# VALEURS DE L'ÉTOILE DE RÉFÉRENCE
#Data_Star(3338,T_proxima_stars,D_proxima_stars,R_proxima_stars,F_proxima_stars) #PROXIMA DU CENTAURE


# VALEURS DE L'ÉTOILE DE RÉFÉRENCE: PROXIMA DU CENTAURE	

p = planetes_names('ProximaCen','b')
#print '\n',W1[p] + '\t' + W2[p] + '\t' + P_l[p]

for i in range(0,len(PROXIMA_T)):
	T_proxima_stars[i] = T_s[p]
	D_proxima_stars[i] = D_s[p]
	R_proxima_stars[i] = R_s[p]
	F_proxima_stars[i] = F_s[p]*1e-26 #Unité SI

#print '\n',W1[p] + '\t' + P_l[p]

figure(figsize=(15, 10), dpi=100)
title(u"Radius vs Temperature " + "\n$\lambda$ ="+ str("%.2e"%L)+"m" ,fontsize = f)
xlabel(u"Effective temperature (°K)",fontsize = f)
ylabel(u"Planetary radius (R_Earth)",fontsize = f)
xscale('log')
#yscale('log')
xlim(xmin=100,xmax=7500)
ylim(ymax=3)
grid(True)



#CODAGES DES COURBES D'ISOCONTRASTES 

#k = C[p]
k = C_lim

PROXIMA_R = pow( k * F_proxima_stars * c**2 * D_proxima_stars**2 / (2*pi*h*nu**3) , 1./2 ) * pow(X,1./2) 

plot(PROXIMA_T,PROXIMA_R/R_Terre,markersize=3,color='red',label='Constrast auround proxima = ' + str("%3.e"%k))

#axhline(y = R_pr[p]/R_Terre)
#axvline(x = T_pr[p])

#plot(T_pr[p],R_pr[p]/R_Terre,'k*',markersize=8)


# VALEURS DE L'ÉTOILE DE RÉFÉRENCE: ROSS 128

p = planetes_names('Ross','128')
#print '\n',W1[p] + '\t' + W2[p] + '\t' + P_l[p]


for i in range(0,len(PROXIMA_T)):
	T_proxima_stars[i] = T_s[p]
	D_proxima_stars[i] = D_s[p]
	R_proxima_stars[i] = R_s[p]
	F_proxima_stars[i] = F_s[p]*1e-26 #Unité SI

k = C_lim

PROXIMA_R = pow( k * F_proxima_stars * c**2 * D_proxima_stars**2 / (pi*h*nu**3) , 1./2 ) * pow(X,1./2) 

#plot(PROXIMA_T,PROXIMA_R/R_Terre,markersize=3,color='red',linestyle='--',label='Constrast auround ROSS 128 = ' + str("%3.e"%k))


# VALEURS DE L'ÉTOILE DE RÉFÉRENCE: Soleil a 10pc 	
for i in range(0,len(PROXIMA_T)):
	T_proxima_stars[i] = 5777
	D_proxima_stars[i] = 10*pc
	R_proxima_stars[i] = 696342e3
	F_proxima_stars[i] = F_solaire(L) #Unité SI

k = C_lim
PROXIMA_R = pow( k * F_proxima_stars * c**2 * D_proxima_stars**2 / (pi*h*nu**3) , 1./2 ) * pow(X,1./2) 
plot(PROXIMA_T,PROXIMA_R/R_Terre,markersize=3,color='orange',linestyle='--',label='Constrast at 10 pc auround The SUN  = ' + str("%3.e"%k))


legend(prop={'size':f})
#show()


###########################################################################
'''CATÉGORIES DE CONTRASTE DES EXOPLANÈTES DANS LE PLAN (T,R)'''
###########################################################################
#Rayon_Temperature(20,C,R)
###########################################################################


###########################################################################
'''VISILITÉS DES EXOPLANÈTES EN FONCTIONS DES HÉMISPHÈRES'''
###########################################################################
#Hemisphere(C,f,**data)
###########################################################################


###########################################################################
'''EXOPLANÈTES VISIBLES DEPUIS PARANAL'''
###########################################################################
#Paranal(C,f,R,**data2)
###########################################################################

###########################################################################

###########################################################################
'''FLUX PLANETES VS DISTANCE A L'ETOILE'''
###########################################################################
#Flux_Distance(D_s,F_p)
###########################################################################
