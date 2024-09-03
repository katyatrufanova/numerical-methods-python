import numpy as np
import timeit
import matplotlib.pyplot as plt

#Algoritmo di sostituzione all'indietro
#====================================
#Parametri di input:
# A : matrice triangolare superiore
# b : vettore dei termini noti
#====================================
#Parametri di output:
#x :  vettore soluzione del sistema

def SostituzioneIndietro(A,b):
    n = len(A)
    x = np.zeros(n)
    if abs(np.prod(np.diag(A))) < 1.0e-14:
        print('Errore: matrice possibilmente singolare!')
    else:
        for i in range(n-1,-1,-1):
            S = 0 
            for j in range(i+1,n):
                S = S + A[i,j]*x[j]
            x[i] = (b[i] - S)/A[i,i]
    return x

def SostituzioneIndietroOttimizzato(A,b):
    n = len(A)
    x = np.zeros(n)
    if abs(np.prod(np.diag(A))) < 1.0e-14:
        print('Errore: matrice possibilmente singolare!')
    else:
        for i in range(n-1,-1,-1):
            S = 0 
            S = S + A[i,i+1:n]*x[i+1:n]
            x[i] = (b[i] - S)/A[i,i]
    return x

nmax = 500
nrange = range(5,nmax,10)
TempoIndietro = np.zeros(len(nrange))
TempoIndietroO = np.zeros(len(nrange))
i = 0
for n in nrange:
    # Costruzione di un problema test
    A = (2*np.random.random((n,n))-1)*10
    xsol = np.ones((n,1))
    L = np.triu(A)
    b = np.dot(L,xsol)
    
    # Calcolo tempo per la sostituzione in avanti
    # Uso della funzione timeit del modulo timeit:
    # restituisce la media del tempo di esecuzione di un comando,
    # calcolata sul numero di iterazioni specificato dal parametro number
    TempoIndietro[i] = timeit.timeit('x = SostituzioneIndietro(L,b)',\
                                   globals=globals(),number=5)

    # Calcolo tempo per la sostituzione in avanti ottimizzata
    TempoIndietroO[i] = timeit.timeit('x = SostituzioneIndietroOttimizzato(L,b)',\
                                    globals=globals(),number=5)
    
    i = i + 1
    
# Rappresentazione grafica del tempo di calcolo
plt.figure(1)
plt.semilogy(nrange,TempoIndietro,'r-',label='Sostituzione all\'indietro')
plt.semilogy(nrange,TempoIndietroO,'g-',label='Versione ottimizzata')
plt.xlabel('Dimensione del problema test')
plt.ylabel('Tempo di calcolo')
plt.legend()
plt.title('Confronto prestazioni')
plt.show()