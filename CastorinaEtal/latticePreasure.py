import numpy as np
import matplotlib.pyplot as mp
import sys

def pLattice(Tarray):

	'''
	[Tarray]  = MeV
	'''

	#parameters
	
	pid = 95/180/np.pi**2
	ct = 3.8706
	an = -8.7704
	bn = 3.9200
	cn = 0.
	dn = 0.3419
	t0 = 0.9761
	ad = -1.2600
	bd = 0.8425
	cd = 0.
	dd = -0.0475

	#Tarray = 130. #MeV

	Tc = 154. #MeV
	t = Tarray / Tc

	f = (pid + an/t + bn/t**2 + cn/t**3 + dn/t**4 ) / (1 + ad/t + bd/t**2 + cd/t**3 + dd/t**4 )
	
	pLattice = Tarray**4 /2. * (1 + np.tanh(ct*(t-t0)) ) * f

	return pLattice

def pHRG(Tarray):
	''' 

	input: 
		[T] = MeV

	return:
		pHRG = pLattice_l * (T / T_l)**4 + g(T), eq (19)
		g(T) = T**4 [ a1 * (T-T_l) + a2 / 3 * (T**3 - T_l**3) + 
				a3 / 4 * (T**4 - T_l**4) + a4 / 10 *(T**10 - T_l**10)]

			a1 = 4.654 GeV**-1
			a2 = -879. GeV**-3
			a3 = 8081. GeV**-4
			a4 = 7.039e6 GeV**-10
			T_l = 130. MeV

	energy density electroweak sector:

	e_ew = 3 * p_ew = g_ew * pi**2 / 30 * T**4, g_ew = 14.45

	p_ew  --> pressure electroweak sector'''

	T = Tarray/1e3 # convert temperature from MeV to GeV

	a1 = 4.654 #GeV**-1
	a2 = -879. #GeV**-3
	a3 = 8081. #GeV**-4
	a4 = -7.039e6 #GeV**-10
	T_l = 0.130 #GeV

	g = T**4 * ( a1 * (T-T_l) + a2 / 3 * (T**3 - T_l**3) + a3 / 4 * (T**4 - T_l**4) + a4 / 10 *(T**10 - T_l**10))

	#evalation pLattice at 130MeV
	#T_lattice = 0.130 #GeV
	pLattice_l = pLattice(130.)

	pressureHRG = pLattice_l * (T / T_l)**4 + g

	return pressureHRG

def energy_ew(Tarray):

	''' 
	energy density e_ew = 3 * p_ew = g_ew * np.pi**2 / 30 * T**4
	'''

	g_ew = 14.45

	return g_ew * np.pi**2 /30 * Tarray**4

def energy_dm(T, energy):
	'''
	energy == energy of non-cold-dark-matter contributions
	[T] = MeV
	'''
	return 2e-10 * (energy)

T_ini = 70.	#MeV
T_end = 400.#MeV
size = 100 #1000 points
step = (T_end - T_ini) / size

T = np.arange(T_ini,T_end, step)

#Lattice (strong?)

#pressureS = pLattice(T)
pressureS = pLattice(T)
pressureS_at130 = pressureS / pLattice(130.)
energyS = 3. / pressureS 
# RGH (electro weak?)

#energyEW = energy_ew(T) 
#pressureEW_at130 = 3. / energy_ew(130.)
#pressureEW = 3. / energyEW /pressureEW_at130

#energyWOcdm = energyS + energyEW

#CDM
#energyCDM = energy_dm(T, energyWOcdm)
#pressureCDM = np.zeros(np.shape(pressureEW)) 

print 'strong E', energyS
#print 'electroweak E', energyEW
#print 'cdm E', energyCDM
print 'strong P', pressureS
#print 'electroweak P', pressureEW
#print 'cdm P', pressureCDM

w_s = pressureS_at130 / energyS
#w_ew = (pressureS + pressureEW) / (energyS + energyEW)
#w_cdm = (pressureS + pressureEW + pressureCDM) / (energyS + energyEW + energyCDM)


mp.plot(T, w_s, 'b-', label = 's')
#mp.plot(T, w_ew, 'r-', label = 's+ew')
#mp.plot(T,w_cdm, 'k-', label = 's+ew+cdm')
mp.show()

