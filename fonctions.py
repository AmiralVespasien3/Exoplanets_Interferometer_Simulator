#!/usr/bin/env python
# coding: utf-8
# programme.py
# variables.py

from pylab import *
from variables import *


#################################################################################
'''TRIAGE croissant dune liste et modification correspondante de la liste associée'''
#################################################################################
def sort2L(A,B):

	i = 0
	s = 0
	while i <= len(A)-2 :

		print i,'\t',s
		if i==300:
			break

		if A[i] > A[i+1]:

			k = A[i]
			A[i] = A[i+1]
			A[i+1] = k

			k = B[i]
			B[i] = B[i+1]
			B[i+1] = k

			i = 0 #REMISE A ZÉRO
			s += 1

		else:

			i += 1

#################################################################################
'''TRIAGE croissant dune liste et modification correspondante des listes associées'''
#################################################################################
def sort2(A,D):

	i = 0
	while i <= len(A)-2 :

		if A[i] > A[i+1]:

			k = A[i]
			A[i] = A[i+1]
			A[i+1] = k

			#MODIFICATION DE LA PLACE DES ÉLÉMENTS DANS TOUTE LES AUTRES LISTES
			for j in range(0,len(D)):
				k = D[j][i]
				D[j][i] = D[j][i+1]
				D[j][i+1] = k			

			i = 0 #REMISE A ZÉRO

		else:

			i += 1

#################################################################################
'''AFFICHAGE AVEC CHIFFRES SIGNIFICATIFS'''
#################################################################################
def signif(x, digit):
    if x == 0:
        return 0
    return round(x, digit - int(math.floor(math.log10(abs(x)))) - 1)

#################################################################################
'''ALGORIHME PRINCIPAL'''
#################################################################################
def Cibles(L,R,C,F_p,F_rp,Time,T_max,B_max):

	# AFFICHAGE DES CIBLES SELECTIONNÉES

	#Mise a 0 du compteur
	s=0 

	fichier_1 = open("Targets_HI5.tbl","w")
	fichier_2 = open("TARGETS_NAMES.tbl","w")

	fichier_1.write("\n")
	fichier_1.write( "Observation wavelength : " + str("%.2e" %L) + " m" )
	fichier_1.write("\n")
	fichier_1.write( "Bond Albedo: " + str("%.2e" %A_B) )
	fichier_1.write("\n")
	fichier_1.write( "The base line is: " + str("%.2e" %B_max) + " m" )
	fichier_1.write("\n")
	fichier_1.write("Flux is the thermal flux of the planet \nR flux is the reflective flux of the planet")
	fichier_1.write("\n")
	fichier_1.write("The contrast is computed for the total flux")

	#############################################################################

	fichier_1.write("\n\n########################################################\n")

	fichier_1.write( "HOSTS STARS\n")

	#############################################################################
	fichier_1.write("\nNORTHERN HEMISPHERE\n")
	#############################################################################
	fichier_1.write( "Parsec" + '\t' + "\t" +'\t' + "Magnitude" + '\t\t' + 'RA' + '\t\t' + 'DEC' + '\t\t\t' + "Star")

	for i in range(0,len(Nr_p)):

		if math.isnan(P_m[i]) is True:
			P_m[i] = 0

		if P_mt[i] != 'V' and P_mt[i] != 'Kepler-band' :
			P_mt[i] = str(0)

		if math.isnan(C[i]) is False and C[i]>C_lim and (0 < P_d[i] <= 90) and S_pa[i]>=R and  Time[i]<T_max:
			fichier_1.write(
				'\n' 
				+ str("%.2e" % (D_s[i]/pc)) 
				+ "\t" 
				+ str(P_mt[i])  
				+ '\t' 
				+ str("%.2e" % P_m[i]) 
				+ '\t' 
				+ str("%.2e"%D[5][i]) 
			 	+ '\t' 
			 	+ str("%.2e"%D[6][i]) 
			 	+ '\t' 
				+ W1[i])
			s += 1	

	#############################################################################
	fichier_1.write("\n\nSOUTHERN HEMISPHERE\n")
	#############################################################################

	fichier_1.write( "Parsec" + '\t' + "\t" +'\t' + "Magnitude" + '\t\t' + 'RA' + '\t\t' + 'DEC' + '\t\t\t' + "Star")

	for i in range(0,len(Nr_p)):

		if math.isnan(P_m[i]) is True:
			P_m[i] = 0

		if P_mt[i] != 'V' and P_mt[i] != 'Kepler-band' :
			P_mt[i] = str(0)

		if math.isnan(C[i]) is False and C[i]>C_lim and (-90 < P_d[i] <= 0) and S_pa[i]>=R and  Time[i]<T_max:
			fichier_1.write(
				'\n' 
				+ str("%.2e" % (D_s[i]/pc)) 
				+ "\t" 
				+ str(P_mt[i])  
				+ '\t' 
				+ str("%.2e" % P_m[i]) 
				+ '\t' 
				+ str("%.2e"%D[5][i]) 
			 	+ '\t' 
			 	+ str("%.2e"%D[6][i]) 
			 	+ '\t' 
				+ W1[i])
			s += 1

	###################################################################################

	fichier_1.write("\n\n########################################################")

	fichier_1.write( "\n\nEXOPLANETS ORBITING HOSTS STARS\n")

	##################################################################################
	fichier_1.write("\nNORTHERN HEMISPHERE \n")
	###################################################################################

	fichier_1.write( "R_J" + '\t\t\t' + "°K" + '\t\t\t' + "mas" + '\t\t\t' + 'FLux(Jy)' + '\t' + 'R Flux(Jy)' + '\t' + 'Contrast' + '\t')

	for i in range(0,len(D[0])):

		if math.isnan(D[3][i]) is True:
			D[3][i] = 0

		if math.isnan(D[4][i]) is True:
			D[4][i] = 0

		if math.isnan(C[i]) is False and C[i]>C_lim and (0 < D[6][i] <= 90) and D[17][i]>=R and Time[i]<T_max:

			fichier_1.write('\n' 
			 + str("%.2e" % (D[3][i]/RJ)) 
			 +'\t'
			 + str("%.2e" % D[4][i]) 
			 + '\t'  
			 + str("%.2e" % D[17][i]) 
			 +'\t'
			 + str("%.2e" % F_p[i]) 
			 + '\t'
			 + str("%.2e" % F_rp[i]) 
			 + '\t'
			 + str("%.2e" % C[i])
			 + '\t' 
			 + str(D[7][i]) 
			 + '\t' 
			 +  D[13][i])

	##################################################################################
	fichier_1.write("\n\nSOUTHERN HEMISPHERE\n")
	##################################################################################

	fichier_1.write( "R_J" + '\t\t\t' + "°K" + '\t\t\t' + "mas" + '\t\t\t' + 'FLux(Jy)'+ '\t' + 'R Flux(Jy)' + '\t' + 'Contrast' + '\t' )

	for i in range(0,len(D[0])):

		if math.isnan(D[3][i]) is True:
			D[3][i] = 0

		if math.isnan(D[4][i]) is True:
			D[4][i] = 0

		if math.isnan(C[i]) is False and C[i]>C_lim and (-90 < D[6][i] <= 0) and D[17][i]>=R and Time[i]<T_max:

			fichier_1.write('\n' 
			 + str("%.2e" % (D[3][i]/RJ)) 
			 +'\t'
			  + str("%.2e" % D[4][i]) 
			 + '\t'  
			 + str("%.2e" % D[17][i]) 
			 +'\t'
			 + str("%.2e" % F_p[i]) 
			 + '\t'
			 + str("%.2e" % F_rp[i])
			 + '\t'
			 + str("%.2e" % C[i])
			 + '\t'
			 + str(D[7][i]) 
			 + '\t' 
			 +  D[13][i])

	###################################################################################

	'''

	fichier_1.write("\n\n\n\n")
	fichier_1.write("\n ADITIONAL VIKING TARGETS \n")

	###################################################################################

	fichier_1.write("\nNORTHERN HEMISPHERE \n")

	fichier_1.write( "R_J" + '\t\t\t' + "mas" + '\t\t\t' + "°K" + '\t\t\t' + 'FLux(Jy)' + '\t' + 'Contrast' + '\t'
		+  'RA' + '\t\t\t' + 'Dec' + '\t\t\t' + 'm_L' + '\t\t\t' + "Hi5_Days" + '\t' + "VLTi_Days" )

	for i in range(0,len(D[0])):

		if math.isnan(D[3][i]) is True:
			D[3][i] = 0

		if math.isnan(D[4][i]) is True:
			D[4][i] = 0

		if math.isnan(C[i]) is False  and (0 < D[6][i] <= 90) and D[17][i]>=R and Time[i]<T_max:

			if (C[i]>4e-6 and  5.5 < m_L[i]< 6.5) or (C[i]>6e-5 and  8.5 < m_L[i]< 9.5) or (C[i]>3e-7 and  2.5 < m_L[i]< 3.5):

				fichier_1.write('\n' 
				 + str("%.2e" % (D[3][i]/RJ)) 
				 +'\t' 
				 + str("%.2e" % D[17][i]) 
				 +'\t'
				 + str("%.2e" % D[4][i]) 
				 + '\t' 
				 + str("%.2e" % F_p[i]) 
				 + '\t'+ str("%.2e" % C[i])
				 + '\t' 
				 + str("%.2e"%D[5][i]) 
				 + '\t' 
				 + str("%.2e"%D[6][i]) 
				 + '\t' 
				 + str("%.2e"%m_L[i])
				 + '\t\t' 
				 + str("%.2e"%Time[i]) 
				 + '\t'
				 + str("%.2e"%(2*Time[i])) 
				 + '\t' 
				 + str(D[7][i]) 
				 + '\t' 
				 +  D[13][i])

	
	fichier_1.write("\nSOUTHERN HEMISPHERE \n")

	fichier_1.write( "R_J" + '\t\t\t' + "mas" + '\t\t\t' + "°K" + '\t\t\t' + 'FLux(Jy)' + '\t' + 'Contrast' + '\t'
		+  'RA' + '\t\t\t' + 'Dec' + '\t\t\t' + 'm_L' + '\t\t\t' + "Hi5_Days" + '\t' + "VLTi_Days" )

	for i in range(0,len(D[0])):

		if math.isnan(D[3][i]) is True:
			D[3][i] = 0

		if math.isnan(D[4][i]) is True:
			D[4][i] = 0

		if math.isnan(C[i]) is False and  5.5 < m_L[i]< 6.5 and (-90 < D[6][i] <= 0) and D[17][i]>=R and Time[i]<T_max:

			if (C[i]>4e-6 and  5.5 < m_L[i]< 6.5) or (C[i]>6e-5 and  8.5 < m_L[i]< 9.5) or (C[i]>3e-7 and  2.5 < m_L[i]< 3.5) :

				fichier_1.write('\n' 
				 + str("%.2e" % (D[3][i]/RJ)) 
				 +'\t' 
				 + str("%.2e" % D[17][i]) 
				 +'\t'
				 + str("%.2e" % D[4][i]) 
				 + '\t' 
				 + str("%.2e" % F_p[i]) 
				 + '\t'+ str("%.2e" % C[i])
				 + '\t' 
				 + str("%.2e"%D[5][i]) 
				 + '\t' 
				 + str("%.2e"%D[6][i]) 
				 + '\t' 
				 + str("%.2f"%m_L[i])
				 + '\t\t'
				 + str("%.2e"%Time[i]) 
				 + '\t'
				 + str("%.2e"%(2*Time[i])) 
				 + '\t' 
				 + str(D[7][i]) 
				 + '\t' 
				 +  D[13][i])

	'''
	###################################################################################

	#fichier_1.write("\n\n\n\n")
	fichier_2.write("TARGETS DATA \n")

	###################################################################################

	#fichier_2.write( "R_J" + '\t\t\t' + "°K" + '\t\t' + "mas" + '\t\t\t' + 'FLux(mJy)' + '\t' + 'Contrast' + '\t' + "Hi5_Hours" + '\t' + 'Names' )
	

	fichier_2.write( 'Names' )

	figure(figsize=(15, 10), dpi=100)
	
	f =20

	title(u"Constrast vs Angular separation \n" + "$\lambda$ = " + str("%.2e"%L) + " m" ,fontsize = f)
	xlabel(u"Angular separartion (mas)",fontsize = f)
	ylabel(u"Contrast",fontsize = f)

	ylim(ymin=1e-10)

	xscale('log')
	yscale('log')
	grid(True)


	print "Don't panic ! It is computing..."

	compteur = 0
	for i in range(0,len(D[0])):

		if math.isnan(D[3][i]) is True:
			D[3][i] = 0

		if math.isnan(D[4][i]) is True:
			D[4][i] = 0

		if math.isnan(C[i]) is False and C[i]>=C_lim and (-90 <= D[6][i] <= 90) and R <= D[17][i] <= OWA and Time[i]<T_max:

			compteur += 1

			fichier_2.write('\n' 
			 #+ str("%.2e" % (D[3][i]/RJ)) 
			 #+'\t'
			 #+ str("%.1i" %	D[4][i]) 
			 #+ '\t\t'  
			 #+ str("%.2e" % D[17][i]) 
			 #+'\t'
			 #+ str("%.2e" % (F_p[i]/(1e-3))) 
			 #+ '\t'
			 #+ str("%.2e" % C[i])
			 #+ '\t'  
			 #+ str("%.2e"%(Time[i]*24)) 
			 #+ '\t'
			 + str(D[13][i])
			 + '\t' 
			 +  D[7][i])
			 #+ '\t\t\t\t\t\t'
			 #+ M_n[i])

			plot(S_pa[i],C[i],'o',label= D[13][i] + '\t'+ D[7][i] )

		else:
			plot(S_pa[i],C[i],'o',color='grey')

	if compteur < 30:
		legend(prop={'size':10})

#################################################################################################################################

	figure(figsize=(15, 10), dpi=100)
	
	f =20

	title(u"Mass vs Temperature "  ,fontsize = f)
	xlabel(u"Temperature(K)",fontsize = f)
	ylabel(u"Mass(Jupiter mass)",fontsize = f)

	#ylim(ymin=1e-10)

	xscale('log')
	yscale('log')
	grid(True)


	print "Don't panic ! It is computing..."

	compteur = 0
	for i in range(0,len(D[0])):

		if math.isnan(D[3][i]) is True:
			D[3][i] = 0

		if math.isnan(D[4][i]) is True:
			D[4][i] = 0

		if math.isnan(C[i]) is False and C[i]>=C_lim and (-90 <= D[6][i] <= 90) and R <= D[17][i] <= OWA and Time[i]<T_max:

			compteur += 1

			fichier_2.write('\n' 
			 #+ str("%.2e" % (D[3][i]/RJ)) 
			 #+'\t'
			 #+ str("%.1i" %	D[4][i]) 
			 #+ '\t\t'  
			 #+ str("%.2e" % D[17][i]) 
			 #+'\t'
			 #+ str("%.2e" % (F_p[i]/(1e-3))) 
			 #+ '\t'
			 #+ str("%.2e" % C[i])
			 #+ '\t'  
			 #+ str("%.2e"%(Time[i]*24)) 
			 #+ '\t'
			 + str(D[13][i])
			 + '\t' 
			 +  D[7][i])
			 #+ '\t\t\t\t\t\t'
			 #+ M_n[i])

			plot(T_pr[i],M_p[i]/MJ,'o',label= D[13][i] + '\t'+ D[7][i] )

		else:
			plot(T_pr[i],M_p[i]/MJ,'o',color='grey')


	fichier_1.close()
	fichier_2.close()

	if compteur <= 30:
		legend(prop={'size':10})
	
	show()	
	
#######################################################################################################################

	print "There is, for now, " + str(len(Nr_p)) + " known exoplanets"
	print "The number of observable exoplanets, without the OWA constraints, with a contrast superior to: " + str("%.2e" % C_lim) + " is : " + str(s) 
	print "PLease, open Targets_Hi5.tbl : It contains data of the observable exoplanets with the OWA constraints"

	#i=3510
	#print '\n'+ str(P_mt[i]) + "\t" + str(data_n[i,0]) +' ' + str(data_n[i,1])

	#print "Affichage TEST"

	#print data_n[i,0],' ',data_n[i,1],'\t', "F_p " + str("%.4g" % F_p[i]) + " Jy", '\t' ,"F_s " + str("%.4g" % F_s[i])+ " Jy", '\t' , " C " + str("%.4g" %C[i])


#################################################################################
'''FLUX STELLAIRE'''
#################################################################################
#Calcul du flux solaire a 10 pc ; Unité SI
def F_solaire (L):
	c = 3.00E8
	nu = c/L
	f = 2*(h*nu**3/c**2)*pow(exp((h*nu)/(k_b*5777)) -1, -1) * (7e8/(10*pc))**2 *pi
	return f

#Calcul du flux solaire a 1 UA ; Unité SI
def F_sol_1UA (L):
	c = 3.00E8
	nu = c/L
	f = 2*(h*nu**3/c**2)*pow(exp((h*nu)/(k_b*5777)) -1, -1) * (7e8/(1*AU))**2 *pi
	return f

###################################################################################
'''FLUX PLANÈTES SYSTÈMES SOLAIRES À 10 pc''' 
###################################################################################
def F_p_sun (L,R,T):
	c = 3.00E8
	nu = c/L
	f = 2*(h*nu**3/c**2)*pow(exp((h*nu)/(k_b*T)) -1, -1) * (R/(10*pc))**2 *pi
	return f

def F_p_sol_1UA (L,R,T):
	c = 3.00E8
	nu = c/L
	f = 2*(h*nu**3/c**2)*pow(exp((h*nu)/(k_b*T)) -1, -1) * (R/(1*AU))**2 *pi
	return f

##################################################################################
'''CONTRASTES FLUX PLANÈTES SYSTÈMES SOLAIRES À 10 pc'''
##################################################################################
def C_sun(L,R,T):
	c = F_p_sun(L,R,T)/F_solaire (L)
	return c

###################################################################################
def F_s	(L,*T):
	c = 3.00E8
	nu = c/L 
	a = np.zeros(len(T))
	for i in range(0,len(T)):
		a[i] = 1E26*2*(h*nu**3/c**2)*pow(exp((h*nu)/(k_b*T_s[i])) -1, -1) * (R_s[i]/D_s[i])**2 *pi
	return a

###################################################################################
def F_s_prox(N,L,T_pr,R_pr,D_pr):
	c = 3.00E8
	nu = c/L 
	a = np.zeros(len(N))
	for i in range(0,len(N)):
		a[i] = 1E26*2*(h*nu**3/c**2)*pow(exp((h*nu)/(k_b*T_pr[i])) -1, -1)* (R_pr[i]/D_pr[i])**2 *pi
	return a

###################################################################################
'''BRILLANCE PLANETAIRE'''
###################################################################################
def B(L,T):
	c = 3.00E8	
	nu = c/L
	return 1E26*2*(h*nu**3/c**2)*pow(exp((h*nu)/(k_b*T)) -1, -1)

###################################################################################
'''FLUX PLANÈTAIRE'''
###################################################################################
def F_pa(L,Nr_p,T_p,R_pr,D_s):

	F_p=np.zeros(len(Nr_p))
	s = 0 #Compteur

	for i in range(0,len(Nr_p)):
		F_p[i] = B(L,T_pr[i])*math.pi*(R_pr[i]/D_s[i])**2
		#print F_p[i]

	return F_p

	#CALCUL DES TROUS DANS UN TABLEAU DE DONNÉES
	#if math.isnan(M_p[i]) is True:
	#	s +=1
	#print F_p

	#print "s = " + str(s)
###################################################################################