import numpy as np
import sys
'''

eta(t) = \int_{0}^{t} dt' / a(t')

si usamos \dot{a} = da/dt, no queda

\eta(a) = \int_{0}^{a} da/ (a \dot{a}) = \int_{0}^{a} da / (a^{2} H(a))

donde: H(a) =  \dot{a} / a = H_{0} * \sqrt{\Omega_{R,0} a^{-4} + \Omega_{M,0} a^{-3} 
											+ \Omega_{K,0} a^{-2} + \Omega_{\Lambda, 0} }

y relacionamos redshift con factor de escala haciendo: 1+z = 1/a

Por lo tanto:

\eta(z) = 1/H0 * \int_{0}^{1/(1+z)} da / \sqrt{\Omega_{R,0} + \Omega_{M,0} a 
											+ \Omega_{K,0} a^{2} + \Omega_{\Lambda, 0} a^{4}}

'''

#defino eta(z) o  eta(a)

def value_sqrt(array_time, H0, Omega, a = True):
	'''
	Inputs:
		Omega = Densidades, diccionario
		array_time = instante del paso del tiempo a partir del tiempo conforme o el redshift
		z = corrimiento al rojo (bool) -se podria implementar-
		a = factor de estacal (bool)
		default: z = True,
				 a = False
	
	Return:

		devuelve el valor de eta(z) o eta(a) para un dado valor de z o a.
	'''	
	sqrt_value = []
	for each in array_time:
		sqrt_value.append(1./np.sqrt(Omega['rad'] + Omega['mat'] * each + Omega['k'] * each**2 + Omega['lambda'] * each**4))

	return sqrt_value

#defino el intervalo de redshift
z = np.linspace(1100,1e-15, 1e5)

#defino valores de los parametros
Omega = {'rad' : 9.24e-5,
		 'mat' : 0.315,
		 'k' : 0.,
		 'lambda' : 0.685}

H0 = 67.3 # km seg-1 Mpc-1
c = 3e5 #km seg-1

#paso a los valores de 'a' para calcular 
a = 1/(1+z)

integral_values = value_sqrt(a, H0, Omega)

eta_value = 1./H0 * np.trapz(integral_values, a)

#mpc2km = 3.086e19 #queda en km
#sec2yr = 1/(3.154e7) 
#print a
print eta_value 


