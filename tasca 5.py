# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 10:49:01 2021

@author: guill
"""

import numpy as np
import matplotlib.pyplot as plt

#%%

punts = 500

def Mandelbrot(n):  # Obtenim el mandelbrot per n iteracions
    x = np.linspace(-2, 1, punts)
    y = np.linspace(-1.5, 1.5, punts)
    cRe, cIm = np.meshgrid(x, y)
    
    pla_c = cRe + 1j * cIm  # c
    z = 0  # Z_0
    VC = 0
    
    for n in range(n):
        mascara = np.abs(z) < 1e10  # Detecció nombres massa grans
        z = mascara * z + (1 - mascara) * 1e10
        z = z ** 2 + pla_c # Z_n+1
        CM = np.abs(z) < 2 # Condició mandelbrot
        VC = VC + CM  # Velocitat de convergència
    return CM, VC



def Julia(n):
    x = np.linspace(-2, 2, punts)
    y = np.linspace(-1.5, 1.5, punts)
    cRe, cIm = np.meshgrid(x, y)
    
    pla_c = cRe + 1j * cIm # c
    z = np.sqrt(pla_c + 0.7269 - 0.1889j) # Z_0
    VC = 0
    
    for n in range(n):
        mascara = np.abs(z) < 1e10  # Detecció nombres massa grans
        z = mascara * z + (1 - mascara) * 1e10
        z = z ** 2 - 0.7269 + 0.1889j # Z_n+1
        CJ = np.abs(z) < 2 # Condició mandelbrot
        VC = VC + CJ  # Velocitat de convergència
    return CJ, VC

#%%

plt.figure(1, figsize=(8,8))

for i in range(1):
    iteracions = int(i ** 1.3) + 1  # petita fórmula per determinar les iteracions en que es calcula
    resultat_mandelbrot = Mandelbrot(iteracions)  # càlcul del mandelbrot
    Punts_Mandelbrot = resultat_mandelbrot[0]  # Triem l'array bidimensional que conté els punts
    
    plt.subplot(3, 3, i+1, title='iteració #' + str(iteracions))
    plt.axis('off')
    plt.imshow(Punts_Mandelbrot, cmap='gray')
    
plt.tight_layout()
plt.show()

#%%

plt.figure(2, figsize=(8,8))

for i in range(9):
    iteracions = int(i ** 1.3) + 1  
    resultat_mandelbrot = Mandelbrot(iteracions) 
    VelConv_Mandelbrot = resultat_mandelbrot[1]  # Triem l'array bidimensional que conté la velocitat de convergència
    
    plt.subplot(3, 3, i+1, title='iteració #' + str(iteracions))
    plt.axis('off')
    plt.imshow(VelConv_Mandelbrot, cmap='jet')

plt.tight_layout()
plt.show()

#%%

plt.figure(3, figsize=(8,8))

for i in range(9):
    iteracions = int(i ** 2.4) + 1  
    resultat_julia = Julia(iteracions)  # càlcul de Júlia
    Punts_Julia = resultat_julia[0]  # Triem l'array bidimensional que conté els punts
    
    plt.subplot(3, 3, i+1, title='iteració #' + str(iteracions))
    plt.axis('off')
    plt.imshow(Punts_Julia, cmap='gray')
    
plt.tight_layout()
plt.show()

#%%

plt.figure(4, figsize=(8,8))

for i in range(9):
    iteracions = int(i ** 2.4) + 1  
    resultat_julia = Julia(iteracions)  
    VelConv_Julia = resultat_julia[1]  # Triem l'array bidimensional que conté la velocitat de convergència
    
    plt.subplot(3, 3, i+1, title='iteració #' + str(iteracions))
    plt.axis('off')
    plt.imshow(VelConv_Julia, cmap='jet')    

plt.tight_layout()
plt.show()


