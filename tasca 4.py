# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 11:20:08 2021

@author: guill
"""

import numpy as np
import scipy.integrate as spi
import scipy.special as special
import matplotlib.pyplot as plt


#%%

# 4.1.1.PÈNDOL SIMPLE

g = 9.81  # gravetat
L = 1  # longitud

def equdif(y,t):
    """
    y[0]=velocitat
    y[1]=angle
    """
    return -g * y[1] / L, y[0]

punts = 100
t = np.linspace(0, 10, punts) # array amb els punts del temps
ci = [0, 0.01*np.pi]

res = spi.odeint(equdif, ci, t) # array amb les solucions de l'eq. dif.

plt.figure() # representació dels resultats

plt.subplot(221)
plt.plot(t, res[:,1])
plt.xlabel('Temps(s)')
plt.ylabel('Angle (rad)')

plt.subplot(223)
plt.plot(t, res[:,0], color='orange')
plt.xlabel('Temps(s)')
plt.ylabel('Velocitat (rad/s)')

plt.show()
#%%

# 4.1.2.PÈNDOL 

def equdif1 (y,t):
    # omega = y[0]
    # theta = y[1]
    return -g * np.sin(y[1]) / L, y[0]

CI_pendol = [0, 0.1*np.pi]

resultat1 = spi.odeint(equdif1, CI_pendol, t) 

plt.figure() # representació dels resultats

plt.subplot(221)
plt.plot(t, resultat1[:,1])
plt.xlabel('Temps(s)')
plt.ylabel('Angle (rad)')

plt.subplot(223)
plt.plot(t, resultat1[:,0], color='orange')
plt.xlabel('Temps(s)')
plt.ylabel('Velocitat (rad/s)')

plt.show()

#%%

# 4.1.3.PÈNDOL FORÇAT

m = 0.5  # massa
b = 0.9 # ctt de fregament

def equdif2 (y,t):
    # omega = y[0]
    # theta = y[1]
    return -g * np.sin(y[1]) / L - b * y[0] / L**2 / m , y[0]

resultat2 = spi.odeint(equdif2, CI_pendol, t)

plt.figure()

plt.subplot(221)
plt.plot(t, resultat2[:,1])
plt.xlabel('Temps(s)')
plt.ylabel('Angle (rad)')

plt.subplot(223)
plt.plot(t, resultat2[:,0], color='orange')
plt.xlabel('Temps(s)')
plt.ylabel('Velocitat (rad/s)')

plt.show()

#%%

# 4.1.4.PÈNDOL ESMORTEÏT I FORÇAT

A = 1.35
m_2 = 1
g_2 = 1
L_2 = 1
b_2 = 0.5
Omega = 0.666

def equdif3 (y,t):
    # omega = y[0]
    # theta = y[1]
    return -g_2 * np.sin(y[1])/L_2 + (-b_2 * y[0] + A * np.cos(Omega * t))/L_2**2/m , y[0]

resultat3 = spi.odeint(equdif3, CI_pendol, t)

plt.figure()

plt.subplot(221)
plt.plot(t, resultat3[:,1])
plt.xlabel('Temps(s)')
plt.ylabel('Angle (rad)')

plt.subplot(223)
plt.plot(t, resultat3[:,0], color='orange')
plt.xlabel('Temps(s)')
plt.ylabel('Velocitat (rad/s)')

plt.show()

#%%

# 4.2.CONDICIONS INICIALS DEFINIDES EN PUNTS ARBITRARIS 
# SCHRÖDINGUER

def Hermite (n,x): # Solució analítica de l'eqüació diferencial
    hermite = 1/np.sqrt((2**n)*np.math.factorial(n)*np.sqrt(np.pi))*np.e**(-x**2/2)*special.eval_hermite(n, x)
    return hermite
    
def equdif4 (y,x): # Eq. dif. a resoldre
    # Psi' = y[0]
    # Psi = y[1]
    return  -(((2*n) + 1) - B*x**2) * y[1], y[0]  # Psi'', Psi'

B = 1
CI_schrod = np.array([[0, 0.751126], [0, -0.531125], [0, 0.459969]])


plt.figure() # representació solució analítica

plt.subplot(221)
x = np.linspace(-5, 5, 1000) # punts de càlcul (posició)

for n in [0, 2, 4]: 
    plt.plot(x,Hermite(n,x))
    plt.xlabel('Posició (x)')
    plt.ylabel('Hermite')

plt.axis([-5, 5, -0.6, 0.8]) # ajust arbitrari dels eixos

plt.subplot(222) # representació solució eqdif de schrödinguer
x = np.linspace(0,5,500) 
# la meitat dels punts de càlcul (posició) ja que ens donen les condicions inicials a x = 0

for n in [0, 2, 4]:
    y_p = np.array(spi.odeint(equdif4, CI_schrod[int(n/2)], x))
    y_n = np.array(spi.odeint(equdif4, CI_schrod[int(n/2)], x[:-1]))
    y = np.append(y_n[::-1, 1],y_p[:, 1])
    plt.plot(np.linspace(-5, 5, 999), y)
    plt.xlabel('Posició (x)')
    plt.ylabel('Schrödinguer')

plt.axis([-5, 5, -0.6, 0.8])
plt.tight_layout() # Marges més estrets
plt.show()

#%%

# 4.3.OSCIL·LADORS ACOBLATS

def equdif5(y, t):
    v1, x1, v2, x2 = y[0], y[1], y[2], y[3] # assignació
    
    # equacions '_d' indica la derivada del valor
    v1_d = (k2 * x2 - (k1 + k2) * x1)/m
    x1_d = v1
    v2_d = (k2 * x1 - (k1 + k2) * x2)/m
    x2_d = v2
    return v1_d, x1_d, v2_d, x2_d

# condicions de contorn
k1, k2 = 10, 0.5
CI_acoblat = [0,1,0,0] 
t_acoblat = np.linspace(0, 40, 1000)
resultat5 = spi.odeint(equdif5, CI_acoblat, t_acoblat)

plt.figure() # representació 

plt.subplot(211)
plt.plot(t_acoblat, resultat5[:,1])
plt.plot(t_acoblat, resultat5[:,3])
plt.xlabel('Temps(s)')
plt.ylabel('Posició objecte')

plt.subplot(212)
plt.plot(t_acoblat, resultat5[:,0])
plt.plot(t_acoblat, resultat5[:,2])
plt.xlabel('Temps(s)')
plt.ylabel('Velocitat objecte')

plt.tight_layout()
plt.show()

#%%

# 4.4.MOVIMENT D'UN COS EN UN SISTEMA GRAVITATORI

def equdif6(y, t):
    vx, x, vy, y = y[0], y[1], y[2], y[3] # assignació
    
    # equacions '_d' indica la derivada del valor
    vx_d = -G*M*x/(x**2+y**2)**(3/2)
    x_d = vx
    vy_d = -G*M*y/(x**2+y**2)**(3/2)
    y_d = vy
    return vx_d, x_d, vy_d, y_d

# Condicions de contorn (Halley-Sol) 
G, M = 6.67e-11, 1.9891e30
au = 1.49598e11
aphelion_H = 35.082 * au
velap_H = 0.869e3 
Periode_H = 3.15576e7 * 75.32

CI_halley = [0, aphelion_H, velap_H, 0] 
t_grav = np.linspace(0, Periode_H, 1000)
resultat6 = spi.odeint(equdif6, CI_halley, t_grav)

plt.subplot(221)
plt.plot(t_grav, resultat6[:,1], color='deepskyblue')
plt.xlabel('Temps (s)')
plt.ylabel('Posició objecte X(m)')

plt.subplot(222)
plt.plot(t_grav, resultat6[:,3], color='royalblue')
plt.xlabel('Temps (s)')
plt.ylabel('Posició objecte Y(m)')

              # mòdul de la velocitat 
velocitat_H = (np.abs(resultat6[:,0]**2 + resultat6[:,2]**2)) ** 0.5

plt.subplot(223)
plt.plot(t_grav, velocitat_H, color='orange')
plt.xlabel('Temps (s)')
plt.ylabel('Velocitat objecte (m/s)')

plt.subplot(224)
plt.plot(resultat6[:, 1], resultat6[:, 3], color='seagreen')
plt.xlabel('x (m)')
plt.ylabel('y (m)')

plt.tight_layout()
plt.show()

# condicions de contorn (Terra-Sol) 
aphelion_T = 152.099e9
velap_T = 29.29e3
Periode_T = 365.256 * 24 * 3600

CI_Terra = [0, aphelion_T, velap_T, 0] 
t_grav = np.linspace(0, Periode_T, 1000)
resultat7 = spi.odeint(equdif6, CI_Terra, t_grav)

plt.subplot(221)
plt.plot(t_grav, resultat7[:,1], color='deepskyblue')
plt.xlabel('Temps (s)')
plt.ylabel('Posició objecte X(m)')

plt.subplot(222)
plt.plot(t_grav, resultat7[:,3], color='royalblue')
plt.xlabel('Temps (s)')
plt.ylabel('Posició objecte Y(m)')

velocitat_T = (np.abs(resultat7[:,0]**2 + resultat7[:, 2]**2))**0.5

plt.subplot(223)
plt.plot(t_grav, velocitat_T, color='orange')
plt.xlabel('Temps (s)')
plt.ylabel('Velocitat objecte (m/s)')

plt.subplot(224)
plt.plot(resultat7[:,1], resultat7[:,3], color='seagreen')
plt.xlabel('x (m)')
plt.ylabel('y (m)')

plt.tight_layout()
plt.show()