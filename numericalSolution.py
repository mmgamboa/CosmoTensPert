import numpy as np
from scipy.integrate import odeint
from scipy.special import j0, j1, y1 # bessel function
import matplotlib.pyplot as plt
import sys

'''

Este archivo se compromete a resolver la ecuacion de onda:
	h'' + 2 * a'/a * h + k^2 * h = 0 , valido para ambos planos de polarizacion.

'''

def fbessel(H, eta, k, mat = False):
    '''
	Cambio de variable. Paso de ecuacion de segundo orden a primer orden para usar integrador odeint().

	H = [1.,0.], condiciones iniciales
	eta = variable de integracion (tiempo conforme)

	Por default calcula para RD (si mat = True -> calcula para MD)
    '''
    h = H[0]
    g = H[1]
    # c is a constant : c = 1 for RD , c = 2 for MD
    c = 1.
    if mat:
    	c = 2.
  	# llamo h' = g
  	# despejo h'' de la ec de onda.
    dhdEta = g
    dgdEta = 1.0 / eta**2 * (-2.* c * eta * g - k**2 * eta**2  * h)

    return [dhdEta, dgdEta]


'''
Existen tres regimenes de la fucion de transferencia: 
	i  )	(tau < tau_dec, k > k_eq), la funcion de transferencia correspondiente a los modos que reingresan al horizonte en RD
	ii )	(tau > tau_dec, k > k_eq), la funcion de transferencia correspondiente a los modos que reingresan al horizonte en RD 
			que luego lo afecta en forma diferente cuando el universo pasa a estar MD
	iii)	(tau, k < k_eq), la funcion de transferencia correspondiente a los modos que reingresan al horizonte en MD;

	i)  --> transfRad
	ii) --> transfRadMat
	iii)--> transfMat
'''
def transfRad(x, k):

	return j0(k*x)

def A(k, eta_eq):

	aux1 = 3./(2*k*eta_eq)
	aux2 = np.cos(2*k*eta_eq) / (2*k*eta_eq)
	aux3 = np.sin(2*k*eta_eq) / (k*eta_eq)**2

	return aux1 - aux2 + aux3

def B(k, eta_eq):

	aux1 = 1. / (k*eta_eq)**2 
	aux2 = np.cos(2*k*eta_eq) / (k*eta_eq)**2
	aux3 = np.sin(2*k*eta_eq) / (2*k*eta_eq)

	return -1 + aux1 - aux2 - aux3

def transfRadMat(x, k, k_eq, eta_eq):
	'''
	x = tiempo conforme
	k = modo (fijo)
	k_eq = tamanio del horizonte al momento de desacople
	eta_eq = tiempo conforme del desacople'''

	a = A(k, eta_eq)
	b = B(k, eta_eq)

	term1 = a * j1(k*x)
	term2 = b * y1(k*x)

	return eta_eq / x * (term1 * term2)

def transfMat(x, k):

	return (3*j1(k*x)/(k*x)) 
	
def turnerEtAll(x,k,hubble):

	y =  k * x / hubble / 370

	return np.sqrt(1 + 1.34*y + 2.5*y**2)

# valor aproximado del parametro de Hubble
hubble = 0.7

# k_eq. Sale de que el tamanio del horizonte al momento de desacople cumple que:
# k_eq * eta_eq = 1, si eta_eq = 0.02 --> k_eq = 1/0.02
eta_eq = 0.02
k_eq = 1 / eta_eq # \sim 50

# No puede integrarse a partir de cero, pues errores numericos. Ponemos un valor inicial que represente ese 'cero'
x0 = 1e-15

''' Condiciones iniciales para resolver la ecuacion
que h/h_prim(0) = 1 es por la normalizacion'''
y0 = 1	
z0 = 0
Y0 = [y0, z0]

# Intervalo temporal de integracion con un millon de puntos equidistantes. Normalizado a eta0 (pues llega hasta 1)
xspan = np.linspace(1e-15, 1.,1e5)
#Normalizacion del tiempo conforme. Asumo eta0 = 1 (corregir luego)
eta0 =1.

# Calculo las funciones de transferencia para los diferentes regimenes...

transfer = []

#
# ================================================= Figura 6 de Watanabe&Komatsu ===============================
# Esta figura representa las soluciones numericas de la ecuacion de onda para las GW primordiales. 
# ============ Modos que quiero conocer
k1 = 20.#*eta_eq
k2 = 200.#*eta_eq
k3 = 1000.#*eta_eq
k_arr = [k1, k2, k3]
# Resuelvo numericamente la ecuacion diferencial

numerical_sol = []
test = [True, False, False] # because k=20 is purely inside the horizon in MD, the other are ambiguos so we take then reenter in RD.
for i, each in enumerate(k_arr):
	numerical_sol.append(odeint(fbessel, Y0, xspan, args=(each,test[i])))

print 'numerical sol shape', np.shape(numerical_sol)
print 'transf shape', np.shape(transfer)

plt.clf()
plt.figure(figsize=(10,8))
plt.xscale('log')
plt.xlim(0.0001,1.)
plt.ylim(-0.2,1.2)
plt.xlabel('$\eta/\eta_{0}$')
plt.ylabel('$h/h^{prim}$')
color = ['b','g','r']
line_type = [':', '--', '-']
for i, each in enumerate(k_arr):
	plt.plot(xspan, numerical_sol[i][:,0], line_type[i], label='numerical sol k={}'.format(each), color = color[i])
	plt.axvline(1./each, ls = line_type[i], color = color[i])

plt.legend(loc = 'upper right', fontsize='medium')
plt.show()
#plt.savefig('figure6-watanabe-numerical-solutions')


#
# ================================================= Figura 7 de Watanabe&Komatsu ===============================
# Esta figura resulta para la comparacion entre las soluciones numericas y analiticas de las funciones de 
# transferencia para dos modos: un modo que claramente entra durante RD (k=1000) al horizonte y otro que lo hace 
# cuando el universo es MD (k=20)
# ============ Modos que quiero conocer


# ==============

#for eta_i in xspan:
#	if k_arr[0]*eta_i < 1.:
#		transfer.append(1.)
#
#	elif k_arr[0]*eta_i > 1.:
#		if eta_i < eta_eq:
#			transfer.append(transfRad(eta_i, k_arr[0]))
#
#		if eta_i > eta_eq:
#			transfer.append(transfRadMat(eta_i, k_arr[0], k_eq, eta_eq))

#for i, k in enumerate(k_arr):
#	aux1 = [] 
#	if k > k_eq*eta_eq:
#		for each in xspan:
#			if each < eta_eq:
#				aux1.append(transfRad(each, k))
#			elif each > eta_eq:
#				aux1.append(transfRadMat(each, k, k_eq, eta_eq))
#
#	if k < k_eq * eta_eq:
#		for each in xspan:
#			aux1.append(transfMat(each, k))
#
#	transfer.append(np.array(aux1))

# ============ Modos que quiero conocer

eta = np.linspace(1e-8, 1, 1e4) #seria el eta sin normalizar a eta0

k1 = 20.#*eta_eq
k2 = 1000.#*eta_eq
k_arr = [k1, k2]
# Resuelvo numericamente la ecuacion diferencial

numerical_sol_2 = []
test = [True, False]
for i, each in enumerate(k_arr):
	numerical_sol_2.append(odeint(fbessel, Y0, eta, args=(each,test[i])))

mat_t = transfMat(eta, k_arr[0])
rad_t = transfRad(eta, k_arr[1])
plot_analytic = [mat_t,rad_t]

plt.clf()
plt.figure(figsize=(10,8))
#plt.xscale('log')
plt.xlim(0.0001,20)
#plt.ylim(-0.2,1.2)
plt.xlabel('$k \eta$')
plt.ylabel('$h/h^{prim}$')
color = ['b','g','r']
line_type = [':', '--', '-']
# ploting numerical solution
for i, each in enumerate(k_arr):
	plt.plot(each*eta, numerical_sol_2[i][:,0]/max(numerical_sol_2[i][:,0]), line_type[1], 
			label='numerical sol k={}'.format(each), color = color[i])
#ploting analytical solution
for i, each in enumerate(k_arr):
	plt.plot(each*eta, plot_analytic[i]/max(plot_analytic[i]), ls = line_type[2], label='anal sol k={}'.format(each), color = color[i])

plt.legend(loc = 'upper right', fontsize='medium')
plt.show()

# ==============
