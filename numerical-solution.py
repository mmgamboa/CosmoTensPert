import numpy as np
from scipy.integrate import odeint
from scipy.special import j0,j1 # bessel function
import matplotlib.pyplot as plt
import sys

def fbessel(H, eta, k):
    
    h = H[0]
    g = H[1]
  
    dhdEta = g
    dgdEta = 1.0 / eta**2 * (-2.* eta * g - k**2 * eta**2  * h)
    return [dhdEta, dgdEta]


def transferfunctionRad(x, k):

	return j0(k*x) / j0(k*0) * ((k*x) / 3 / j1(k*x))

def transferfunctionMat(x, k):

	return (3*j1(k*x)/(k*x)) / (1) * ((k*x) / 3 / j1(k*x))
#	return (3*j1(k*x)/(k*x)) / (j1(k*3e-4)/(k*x)) * ((k*x) / 3 / j1(k*x))

def turnerEtAll(x,k,hubble):

	y =  k * x / hubble / 370

	return np.sqrt(1 + 1.34*y + 2.5*y**2)

hubble = 0.7

x0 = 1e-15
y0 = 1
z0 = 0
Y0 = [y0, z0]
 
xspan = np.linspace(1e-15, 1.,1e6)

#Asumo eta0 = 1 (corregir luego)
eta0 =1.
#eta0 = 1.96
#k = 1000./eta0

# calculo la funcion de transferencia para la radiacion y materia
transf = []
x_rad = []
x_mat = []
for each in xspan:
	if each < 0.02:
		x_rad.append(each)
	elif each >= 0.02:
		x_mat.append(each)

x_rad = np.array(x_rad)
x_mat = np.array(x_mat)
#print x_rad[0], x_rad[-1], x_mat[0],x_mat[-1]
k=10.
k2=1000.
# Resuelvo numericamente la ecuacion diferencial
sol = odeint(fbessel, Y0, xspan, args=(k,))
#Y0=[sol[-1][0],z0]
sol2 = odeint(fbessel, Y0, xspan, args=(k2,))
#sol_m = odeint(fbessel, Y0, xspan, args=(k,))
#sol2_m = odeint(fbessel, Y0, xspan, args=(k2,))


transf_rad = transferfunctionRad(xspan, k)
transf_mat = transferfunctionMat(xspan, k)

# calculo la funcion ajustada por Turner, White and Lindsey (1993)
mod = turnerEtAll(xspan,k,hubble)

plt.clf()
plt.figure(figsize=(10,8))
plt.xscale('log')
plt.xlim(0.0001,1.)
plt.ylim(-0.2,1.2)
plt.xlabel('$\eta$')
plt.ylabel('h')
plt.plot(xspan, sol[:,0], '-', label='numerical sol k={}'.format(k), color = 'b')
plt.plot(xspan, sol2[:,0], '-', label='numerical sol k={}'.format(k2), color = 'k')
#plt.plot(x_rad, transf_rad, '-', label = 'Transfer function (RD)', color = 'g')
#plt.plot(x_mat, transf_mat, '-', label = 'Transfer function (MD)', color = 'k')
plt.axvline(x = 0.02, label = '$\eta_{eq}$=0.02 ', color = 'r', ls = '--')
#plt.plot(xspan, mod, '-', label = 'Turner et all', color = 'y')

plt.legend(loc = 'upper right', fontsize='medium')
plt.show()
#plt.savefig('')