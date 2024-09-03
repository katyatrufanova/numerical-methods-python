import numpy as np

#Algoritmo di sostituzione in avanti
#====================================
#Parametri di input:
# A : matrice triangolare inferiore
# b : vettore dei termini noti
#====================================
#Parametri di output:
#x :  vettore soluzione del sistema

def SostituzioneAvanti(A,b):
    n = len(A)
    x = np.zeros(n)
    if abs(np.prod(np.diag(A))) < 1.0e-14:
        print('Errore: matrice possibilmente singolare')
    else:
        for i in range(0,n):
            S = 0
            for j in range(0,i+1):
                S = S + A[i,j]*x[j]
            x[i] = (b[i] - S)/A[i,i]
    return x

import matplotlib.pyplot as plt

nmax = 1000
nrange = range(5,nmax,10)
Errore = np.zeros(len(nrange))
i = 0
for n in nrange:
    # Costruzione di un problema test
    A = (2*np.random.random((n,n))-1)*10
    xsol = np.ones((n,1))
    L = np.tril(A)
    b = np.dot(L,xsol)
    x = SostituzioneAvanti(L,b)
    errRel = np.linalg.norm(xsol-x)/np.linalg.norm(xsol)
    Errore[i] = errRel
    i = i + 1
    
# Rappresentazione grafica dell'errore
plt.figure(1)
plt.semilogy(nrange,Errore,'r-')
plt.xlabel('Dimensione del problema test')
plt.ylabel('Errore')
plt.title('Errore al variare delle dimensioni')
plt.show()

def SostituzioneAvantiOttimizzato(A,b):
    n = len(A)
    x = np.zeros(n)
    if abs(np.prod(np.diag(A))) < 1.0e-14:
        print('Errore: matrice possibilmente singolare')
    else:
        for i in range(0,n):
            S = 0
            S = S + A[i,0:i+1]*x[0:i+1]
            x[i] = (b[i] - S)/A[i,i]
    return x

