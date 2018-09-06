#!/usr/bin/env python
# coding: utf-8
# programme.py

from pylab import *
import mr_forecast as mr




###########################################################################################
#CONSTANTES FONDAMENTALES
##########################################################################################

#Unités du système international
h = 6.63E-34
c = 3.00E8
k_b = 1.38E-23
pc = 3.08E16

M_so = 1.988E30 #Masse Solaire (kg)
R_so = 6.96E8 # Rayon solaire
R_Terre = 6378.137E3 # mètres
AU = 1.50E11
RJ=7.15E7
MJ= 1.90E27
M_Terre = 5.97e24 # kg
T_Terre = 254.3 #Kelvin
TJ = 110.0 #Kelvin
rho_J = 1326 # kg/m³ Masse volumique globale de Jupiter
rho_Terre = 5515 # kg/m³ Masse volumique globale de la Terre

m2r = 4.85E-9
A_B = 0.3 #Albedo de Bond Jovien

###########################################################################################
#PARAMÈTRES INTERFEROMÉTRIQUES
##########################################################################################

###############################################################################
#CALCUL RESOLUTION ANGULAIRE
###############################################################################

print '\n'
L = float(raw_input("Enter the wavelength of observation in meters (ex: 3.8e-6): "))
print '\n'
C_lim = float(raw_input("Enter the inferior contrast limit (ex: 1e-5): "))
print '\n'
OWA = float(raw_input("Enter the Outer Working Angle in mas (ex: 250) or ignore it by typing '0' : "))
if OWA==0:
	OWA = 1e5
print '\n'
choice = float(raw_input('If you choose the VLTI: Press 1 \n If you choose the LIFE interferometer : Press 2  \n If you want to have free choice: Press 3 \n '))
if choice==1:
	B_max = 200 #Longeur de base max du VLT (m)
elif choice==2:
	B_max = 500
	print "You have chosen 'Life' the contrast limit is reset to 1e-7 and the base line is set to 500 meters "
	C_lim = 1e-7
elif choice==3:
	R = float(raw_input("Enter your angular resolution in mas (ex: 2): "))
	B_max = (2*(m2r*R)/L)**(-1)


R = L/(2*B_max)*(1/m2r) #mas
print "\n \n"+"Inner Working Angle = " + str("%.3g" % R) + " mas"	

nu = c/L #Fréquence d'observation

#Bande L
#L = 3.8E-6
#S = 6E-3/(2*sqrt(24)) # Sensibilité Hi5 en 24h (1 journée)
S= 1e-42
#print "\n\n Sensibilité = " + str(S)

#Bande 5.6µm
#L = 5.6E-6
#S = 0.16E-6/(sqrt((24/9.722))) # Sensibilité Hi5 en 24h (1 journée)
#S= 1e-42
#print "\n\n Sensibilité = " + str(S)


#Bande M
# M band 4.82 0.35 6
#L = 4.82E-6
#S = 6E-3/(2*sqrt(24)) # Sensibilité Hi5 en 24h (1 journée)
#print "S = " + str(S)

#Bande N
#B10.7 10,65 1.37 4
#L = 10E-6 
#S = 4E-3/(2*sqrt(24)) # Sensibilité Hi5 en 24h (1 journée)
#S= 1e-42


#Bande N
#L = 10E-6 
#S = 0.54E-6/(sqrt((24/9.722))) # Sensibilité Hi5 en 24h (1 journée)
#S= 1e-42

#Bande N
#L = 15E-6 
#S = 1.39E-6/(sqrt((24/9.722))) # Sensibilité Hi5 en 24h (1 journée)
#S= 1e-42


'''Sensibilité'''
#print "\n S = ",'%.2e'%S


T_max = 1  # Nombres de jours max acceptables pour intégration du VLT
#OWA   = 250 # Outer working angle in mas
#B_max = 8.0 #Diamètre du Télescope Australien (m)
Trans = 0.5 # Facteur de Transmission des optiques

#############################################################################################
#EXTRACTIONS DES DONNÉES DE LA NASA
##############################################################################################

data = genfromtxt("./data/planets.tbl",dtype= float,skip_header=14,delimiter=None)
data_s = genfromtxt("./data/Etoile_R.tbl",dtype= float,skip_header=13,delimiter=None)

data_m = genfromtxt("./data/Etoile_mag.tbl",dtype= float,skip_header=12,delimiter=None,usecols=(0,2))
data_mt = genfromtxt("./data/Etoile_mag.tbl",dtype= str,skip_header=12,delimiter=None,usecols=(1,3))
data_mL = genfromtxt("./data/mag_3.4um.tbl",dtype= float,skip_header=10,delimiter=None)

data_n = genfromtxt("./data/planets_name.tab",dtype= str,skip_header=1,delimiter=None)
data_d = genfromtxt("./data/planets_dec.tbl",dtype= float,skip_header=13,delimiter=None,usecols=(2,4))

data_rho = genfromtxt("./data/Masse_Volumique.tbl",dtype= float,skip_header=11,delimiter=None)
data_SM = genfromtxt("./data/Stellar_Mass.tbl",dtype= float,skip_header=10,delimiter=None)

methods = genfromtxt("./data/Methods.tbl",dtype = str, skip_header=11,delimiter=None)
Amas1    = genfromtxt("./data/Star_Table_Amas.dat",dtype = str, skip_header=1,delimiter=None,usecols=(0,1,2,3))
Amas2    = genfromtxt("./data/Star_Table_Amas.dat",dtype = float, skip_header=1,delimiter=None)

Spec_Data1 = genfromtxt("./data/Spectral_Type_Data.dat",dtype = str, skip_header=2,delimiter=None)
Spec_Data2 = genfromtxt("./data/Spectral_Type_Data.dat",dtype = float, skip_header=2,delimiter=None) 

#######################################################################
# CHARGEMENT DES DONNÉES SIMBADS AVEC LES MAGNITUDES EN BANDE K
# Toutes les données ne sont pas connues pour la bande K: 238 sur 262
#######################################################################
Simbad_data = genfromtxt("./outputs/Tagets_Simbad_K_band.txt",dtype = str)
Simbad_data2 = genfromtxt("./outputs/Tagets_Simbad_K_band.txt",dtype = float,usecols=(0,2))
#######################################################################
#DATA UPLOADING FROM VIKING AND BOULET TARGETS

def screen_list(a,b,s):

	#Screen the names of targets for Viking and Boulet

	n0_B = s[:,a];n1_B = s[:,b] ;list(n0_B);list(n1_B);n_B = []

	
	for i in range(0,len(n0_B)):
		n_B.append(n0_B[i]+n1_B[i])

	'''
	for i in range(0,len(n0_B)):
		print n_B[i]
	'''

	return n_B

def Compare_list(a,b):
	c = [] #Empty list to contain common targets between Viking et Boulet
	d = [] #Empty list to contain different targets between Viking et Boulet
	# Sorting algorithm
	for i in range(0,len(a)):
		for j in range(0,len(b)):
			if a[i] == b[j]:
				c.append(a[i])
	return c

def Extraction_list(a,b):
	for i in range(0,len(a)):
		for j in range(0,len(b)):
			if a[i] == b[j]:
				b[j] = 0

	c=[]
	for i in range(0,len(b)):
		if b[i] != 0:
			c.append(b[i])

	return c

Viking = genfromtxt("./data/Viking_targets",dtype = str, skip_header=0,delimiter=None)
Boulet = genfromtxt("./data/Boulet_targets",dtype = str, skip_header=0,delimiter=None)

n_V = screen_list(0,1,Viking)
#print '\t'
n_B = screen_list(10,9,Boulet)

V1 = Compare_list(n_V,n_B)
#print V1

# List of targets from Viking not in BOULET list
e_list = Extraction_list(V1,n_V)
#print '\n'
#print e_list

#print V1,'\n',V2

#######################################################################
'''CREATION DES ARRAYS'''
#######################################################################

Nr_p = data[:,0] #Numéro référence des planètes		
O_p = data[:,1]*AU #Orbite planetes
M_p = data[:,2]*MJ #Masse de la planète (kg)

R_p = data[:,3]*RJ #Rayon de la planète (m)
T_p = data[:,5]    #Température planétaire en K
rho  = data_rho[:,1]*1e3 # Masse volumique de la planète en kilogrammes par mètres cubes 
 

P_ad  = data_d[:,0] # Ascension droite des planètes en degrés décimaux
P_d  = data_d[:,1] #Déclinaison des planètes en degrés décimaux

P_l = data_mt[:,0]  # Lettre de la planète
P_mt = data_mt[:,1] # Bande Magnitude apparente de de l'étoile hôte
P_m = data_m[:,1] #Magnitude apparente de l'étoile hôte
m_L = data_mL[:,1]#Magnitude apparente en Bande Wise 3.4µm (la plus proche de la bande L)
M_S = data_SM[:,1]*M_so #Masse des étoiles (kg)


D_s = data_s[:,1]*pc #Distance stellaire en mètres
T_s = data_s[:,2] 	 #Température effective stellaire en K
R_s = data_s[:,3]*R_so #Rayon stellaire en mètres

W1 = data_n #Première partie du nom de la planète
M_n = methods[:,2] # Méthodes de détection pour chaque découvertes
#print M_n

S_p = data[:,4]*m2r*D_s  #Separation étoile-hôte exoplanète en mètres
S_pua = data[:,4]*m2r*D_s/AU # Separation étoile hôte exoplanète en UA
S_pa = data[:,4] #Separation étoile hôte exoplanète en mas

Names1_A = Amas1[:,0]
Names2_A = Amas1[:,1]
Simb_names = Simbad_data[:,1]
Simb_K = Simbad_data2[:,1]

Kmag_A = zeros(len(Names1_A))

Names_A = []
for i in range(0,len(Names1_A)):
	Names_A.append(Names1_A[i] + Names2_A[i])


for i in range(0,len(Names_A)): 
	for j in range(0,len(Simb_names)):
		#print  Names_A[i],Simb_names[j]
		if Names_A[i] == Simb_names[j]:
			Kmag_A[i] = Simb_K[j]


'''
print 'VERIFICATION'
for i in range(0,len(Names_A)):
	print Names_A[i],Kmag_A[i]

'''


Group_A  = Amas1[:,2]
SType_A  = Amas1[:,3]

Vmag_A   = Amas2[:,4]
VJmag_A  = Amas2[:,6]
Dpc_A    = Amas2[:,8]
RA_A     = Amas2[:,12]
DEC_A    = Amas2[:,13]

SType_DAT = Spec_Data1[:,0]
Mass_DAT  = Spec_Data2[:,1]

Jmag_A = zeros(len(VJmag_A))
AGE_A  = zeros(len(VJmag_A))
Mass_A = zeros(len(VJmag_A))
SType_DAT_N = zeros(len(VJmag_A))


##########################################################################
# PREDICTION DE LA MASSE MINIMUM DANS LES ASSOCIATIONS STELLAIRES PROCHES
##########################################################################

for i in range(0,len(Jmag_A)):
	Jmag_A[i] = float(Vmag_A[i]) - float(VJmag_A[i])

	if Group_A[i] == 'AB':
		AGE_A[i] = 149e-3 # GYR
	elif Group_A[i] == 'Ar':
		AGE_A[i] = 60e-3 # GYR
	elif Group_A[i] == 'BP':
		AGE_A[i] = 24e-3 # GYR
	elif Group_A[i] == 'Ca':
		AGE_A[i] = 45e-3 # GYR
	elif Group_A[i] == 'Co':
		AGE_A[i] = 42e-3 # GYR
	elif Group_A[i] == 'EC':
		AGE_A[i] = 11e-3 # GYR
	elif Group_A[i] == 'TH':
		AGE_A[i] = 45e-3 # GYR
	elif Group_A[i] == 'TW':
		AGE_A[i] = 10e-3 # GYR
	elif Group_A[i] == '32':
		AGE_A[i] = 22e-3 # GYR

#REMPLISSAGE DES MASSES AVEC LA TABLE DE RÉFÉRENCE
	for j in range(0,len(SType_DAT)):
		if SType_A[i][0]+SType_A[i][1]+SType_A[i][2] == SType_DAT[j]:
			Mass_A[i] = Mass_DAT[j] # Solar Mass
			break
		else:
			if SType_A[i] == SType_DAT[j]:
				Mass_A[i] = Mass_DAT[j] # Solar Mass
				break
			else:
				pass
			

Table_S = zeros(9)
Table_A = ['AB','Ar','BP','Ca','Co','EC','TH','TW','32']

for i in range(0,len(Names1_A)):
	for j in range(0,len(Table_A)):
		if Group_A[i] == Table_A[j]:
			Table_S[j] +=1

'''
ns = zeros(9)
print 'Group','\t','Ty','\t','M(sol)','\t','Pc','\t','Gyr','\t','Vmag','\t','Kmag','\t','N1','N2'
for i in range(0,len(Names1_A)):
	for j in range(len(Table_A)):
		if Group_A[i] == Table_A[j] and Jmag_A[i]<30 and (SType_A[i][0]=='F' or SType_A[i][0]=='G' or SType_A[i][0]=='K' or SType_A[i][0]=='M') and ( -70.0 <= DEC_A[i] <= +15) and (Mass_A[i] != 0) :
			print Group_A[i],'\t',SType_A[i],'\t',Mass_A[i],'\t',Dpc_A[i],'\t',AGE_A[i],'\t',Vmag_A[i],'\t',Jmag_A[i],'\t',Names1_A[i],Names2_A[i]
			print Names_A[i]
			ns[j] +=1
#print 's = ',s
'''

#SCREEN SHOT OF THE NUMBER OF STARS IN THE MOOVING GROUPS
'''
for i in range(0,len(Table_S)):
	print "Le nombre d'étoiles dans l'association stellaire",Table_A[i],"est:",Table_S[i]		
	print "Le nombre d'étoiles totales dans l'association stellaire",Table_A[i],"ayant une magnitude inférieure à 10, en bande J est : ",ns[i]
	print '\n'
'''

#######################################################################
# VARIABLE DU CALCUL D'ISOCONTRASTE
#######################################################################

PROXIMA_T = arange(100,7500,10)

X = zeros(len(PROXIMA_T))


T_proxima_stars = zeros(len(PROXIMA_T))
D_proxima_stars = zeros(len(PROXIMA_T))
R_proxima_stars = zeros(len(PROXIMA_T))
N_proxima_stars = zeros(len(PROXIMA_T))
F_proxima_stars = zeros(len(PROXIMA_T))


X = exp( (h*nu)/(k_b*PROXIMA_T) ) - 1

###############################################################################
#MISE A JOUR DES DONNÉES MANQUANTES
###############################################################################

def planetes_names(n1,n2):
	p = 0
	for i in range(0,len(W1)):
		if W1[i]==n1 and P_l[i]==n2:
			p = i
	return p

#Beta pictoris 
p = planetes_names('betPic','b')
R_p[p] = 1.5*RJ

#LkCa 15 c
#Pas de données

#LkCa 15 b
#Pas de données

#Kepler-417 c 
p = planetes_names('Kepler-417','c')
R_p[p] = 0.240*RJ

#Kepler-416 c 
#Données douteuses de Mortan et al.2016

#Kepler-415 c 
p = planetes_names('Kepler-415','c')
R_p[p] = 0.23*RJ

#Kepler-37 e 
#Pas de données

#Kepler-1653 b
p = planetes_names('Kepler-1653','b')
R_p[p] = 0.194*RJ

#KIC-10001893 d,b,c
#Contreversé!

#HR 8799 e
p = 755 - 1
R_p[p] = 1*RJ
#print "\n" + W1[p]+" "+W2[p]

#HD 189733 b
p = 455 - 1
R_p[p] = 1.138*RJ

#GJ 186 A b
p = planetes_names('Kepler-1653','b')
#print W1[p],P_l[p]

#2MASS J12073346-3932539 b (Planète de Gaël Chauvin)
i = 13-1
T_s[i]=2550
R_s[i]=0.25*R_so
T_p[i]=1150
R_p[i] = 1.5*RJ


#HD102195 b
p = planetes_names('HD102195','b')
M_p[p] = 0.45*MJ
S_p[p] = 0.049*AU
S_pa[p] = S_p[p]/D_s[p]/m2r 


#HD 285507 b
p = planetes_names('HD285507','b')
M_p[p] = 0.917*MJ 
S_p[p] = 0.060*AU
S_pa[p] = S_p[p]/D_s[p]/m2r


#HD 27894 b
p = planetes_names('HD27894','b')
T_s[p] = 4920

#HD 16417 b
p = planetes_names('HD16417','b')
T_s[p] = 5936


#HD 160691 c


#################################################################################
#RAYON DE LA PLANÈTE RECALCULÉ
##################################################################################

#Rmedian, Rplus, Rminus = mr.Mstat2R(mean=1.0, std=0.1, unit='Earth', sample_size=100, classify='Yes')
#print 'R = %.2f (+ %.2f - %.2f) REarth' % (Rmedian, Rplus, Rminus)

R_pr = np.zeros(len(Nr_p))
for i in range(0,len(Nr_p)):

	#Test si le rayon est connu
	if math.isnan(R_p[i]) is False:
		R_pr[i] = R_p[i]

	#Rayon inconnu : Recalcul du rayon 
	else:

	
		if W1[i].find('Kepler') == 0 and (M_p[i]/M_Terre)<= 16 :
			#Relation de Wolgang
			R_pr[i] = pow( M_p[i]/(2.7*M_Terre) , 10./13)*R_Terre 
		

		elif M_p[i]<7*M_Terre:

			if math.isnan(rho[i]) is True:
				R_pr[i] = pow( (3./(4*pi) ) * (M_p[i]/rho_Terre),1./3)
			else:
				R_pr[i] = pow( (3./(4*pi) ) * (M_p[i]/rho[i]),1./3)

		else :
			if math.isnan(rho[i]) is True:
				R_pr[i] = pow( (3./(4*pi) ) * (M_p[i]/rho_J),1./3)
			else:
				R_pr[i] = pow( (3./(4*pi) ) * (M_p[i]/rho[i]),1./3)



####################################################################################
# RECALCUL DE LA TEMPERATURE EFFECTIVE PLANETAIRE
####################################################################################
T_pr  = np.zeros(len(Nr_p))
for i in range(0,len(Nr_p)):

	if math.isnan(T_p[i]) is False:
		T_pr[i] = T_p[i]

	else:
		T_pr[i] = pow(  ( R_s[i]**2 * (1-A_B) ) / (4*S_p[i]**2)  ,   1./4  ) * T_s[i]

######################################################################
# MASSES MANQUANTES##################################
#Mmedian, Mplus, Mminus = mr.Rstat2M(mean=0.1, std=0.01, unit='Jupiter', sample_size=100, grid_size=1e3, classify='Yes')
#print 'M = %.3f (+ %.3f - %.3f) MEarth' % (Mmedian, Mplus, Mminus)
#########################################################################

##########################################################################################################################################



#GJ 86 A b
'''
p = planetes_names('GJ86','b') 
print W1[p],P_l[p]
print M_p[p]/MJ
#M_p[p] = 4.01*MJ
'''

#BD+202457 b
'''
p = planetes_names('BD+202457','b')
print W1[p],P_l[p],'\t',p
print M_p[p]/MJ
'''

#HD110014 c
'''
p = planetes_names('HD110014','b')
print W1[p],P_l[p],'\t',p
print M_p[p]/MJ
'''

#11 Com b
'''
p = planetes_names('11Com','b')
print W1[p],P_l[p],'\t',p
print M_p[p]/MJ
'''

##############################################################################################################################

#LIBRAIRIES
data = {"S_pa_N": [],"C_N": [],"S_pa_S": [],"C_S": []}
data2 = {"S_pa_S": [],"C_S": [],"S_pa_P": [],"C_P": []}

D = {}
D[0]      = Nr_p
D[1]      = O_p
D[2]      = M_p
D[3]      = R_pr
D[4]      = T_pr
D[5]      = P_ad
D[6]      = P_d
D[7]      = P_l
D[8]      = P_mt
D[9]      = P_m
D[10]     = D_s
D[11]     = T_s
D[12]     = R_s
D[13]     = W1
D[15]     = S_p
D[16]     = S_pua
D[17]     = S_pa 