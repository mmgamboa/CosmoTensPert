'''
El objetivo de este programa es reproducir las imagenes correspondientes
a la densidad de onda gravitacional primordial de Watanabe&Komatsu (2006).
(https://arxiv.org/abs/astro-ph/0604176)
'''

import numpy as np
import matplotlib.pyplot as mp
import sys
from scipy.integrate import odeint
from scipy.special import j0,j1 # bessel function

def h2a2_tau(k,a):
	h = 0.7
	Omega_rh2 = 4.15e-5
	Omega_m = 1 - Omega_rh2/(h**2) 
	func = 1e4 * (Omega_rh2 / (a**2) + h**2 / a - Omega_rh2 / a)

	return func

def init_power():
	primordial_value = 1.
	return primordial_value**2

def T_prim(tau,k):
	'''La derivada de T es:
	x = tau * k
	dT/dtau[x] = -(k/x^n) * [ A j_{n+1}(x) + B y_{n+1}(x)]
	'''
