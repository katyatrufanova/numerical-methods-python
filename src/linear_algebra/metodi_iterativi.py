import time
import numpy as np
from sostituzione_indietro import SostituzioneIndietro

#Metodo iterativo di Jacobi
#================================
#Parametri di input:
#A :    matrice dei coefficienti
#b:     vettore dei termini noti
#x0 :   approssimazione iniziale
#toll : tolleranza
#Kmax : numero massimo di iterazioni
#================================
#Parametri di output:
#x1 :   vettore soluzione del sistema
#k : numero dell'iterazione a cui è stata trovata la soluzione
def Jacobi(A,b,x0,toll,Kmax):
    n = len(A)
    x1 = np.copy(x0)
    k = 0       #contatore iterazioni
    stop = 0   #raggiungimento criteri di arresto
    while (not(stop) and (k <= Kmax)):
        for i in range(0,n):
            S = 0
            for j in range(0,n):
                if (j != i):
                    S = S + A[i,j] * x0[j]
            x1[i] = (b[i] - S) / A[i,i]
        r1 = b -(np.dot(A,x1))
        k = k + 1
        stop = (np.linalg.norm(r1) <= toll * np.linalg.norm(b)) and \
                (np.linalg.norm(x1 - x0) / np.linalg.norm(x1) <= toll)
        x0 = x1
    if (not(stop)):
        print('Errore: procedimento non raggiunge tolleranza toll \
              in Kmax iterazioni')
    else:
        return x1, k
    
#Metodo iterativo di Gauss-Seidel
#================================
#Parametri di input:
#A :    matrice dei coefficienti
#b:     vettore dei termini noti
#x0 :   approssimazione iniziale
#toll : tolleranza
#Kmax : numero massimo di iterazioni
#================================
#Parametri di output:
#x1 :   vettore soluzione del sistema
#k : numero dell'iterazione a cui è stata trovata la soluzione
def GaussSeidel(A,b,x0,toll,Kmax):
    n = len(A)
    x1 = np.copy(x0)
    k = 0       #contatore iterazioni
    stop = 0    #raggiungimento criteri di arresto
    while (not(stop) and (k <= Kmax)):
        for i in range(0,n):
             S = 0
             for j in range(0,i):
                 S = S + A[i,j] * x0[j]
             for j in range(i+1,n): 
                 S = S + A[i,j] * x1[j] 
             x1[i] = (b[i] - S) / A[i,i]
        r1 = b -(np.dot(A,x1))
        k = k + 1
        stop = (np.linalg.norm(r1) <= toll * np.linalg.norm(b)) and \
                (np.linalg.norm(x1 - x0) / np.linalg.norm(x1) <= toll)
        x0 = x1
    if (not(stop)):
        print('Errore: procedimento non raggiunge tolleranza toll \
              in Kmax iterazioni')
    else:
        return x1, k
    
# Confronto prestazioni metodi iterativi

def EliminazioneGauss(A,b):
    A = np.copy(A) 
    b = np.copy(b)
    n = len(A)
    for j in range(n-1):
        for i in range(j+1,n):
            m = A[i,j]/A[j,j]
            A[i,j] = 0
            for k in range(j+1,n):
                A[i,k] = A[i,k] - m*A[j,k]
            b[i] = b[i] - m*b[j]
    return A, b

def SoluzioneConMEG(A,b):
    U,c = EliminazioneGauss(A,b)
    x = SostituzioneIndietro(U,c)
    return x
    
# Costruzione di un problema test
n = 1000 # Dimensione della matrice
e = np.ones(n)
e1 = np.ones(n-1)
# Costruzione matrice a diagonale predominanza in senso stretto
A = -4*np.diag(e,0) + np.diag(e1,1) + np.diag(e1,-1) 
xsol = np.ones((n,1))
b = np.dot(A,xsol)

x0 = np.zeros((n,1)) # Approssimazione iniziale
toll = 1.0e-6
Kmax = 1000

# Calcolo dei risultati, confronto con la soluzione test
# e visualizzazione dei risulatati del test
# Uso della funzione numpy.allclose
# numpy.allclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False)
# Returns True if the two arrays are equal within the given tolerance; 
# False otherwise.
xJacobi, k = Jacobi(A,b,x0,toll,Kmax)
xGaussSeidel, k =  GaussSeidel(A,b,x0,toll,Kmax)
print(np.allclose(xsol, xJacobi))
print(np.allclose(xsol, xGaussSeidel))

nmax = 100
nrange = range(5,nmax,10)
toll = 1.0e-6
Kmax = 1000
TempoJacobi = np.zeros(len(nrange))
TempoGaussSeidel = np.zeros(len(nrange))
TempoMEG = np.zeros(len(nrange))
IterazioniJacobi = np.zeros(len(nrange))
IterazioniGaussSeidel = np.zeros(len(nrange))
i = 0
for n in nrange:
    # Costruzione di un problema test
    e = np.ones(n)
    e1 = np.ones(n-1)
    A = -4*np.diag(e,0) + np.diag(e1,1) + np.diag(e1,-1) 
    xsol = np.ones((n,1))
    b = np.dot(A,xsol)
    
    x0 = np.zeros((n,1)) # Approssimazione iniziale
    
    # Calcolo tempo e iterazioni Jacobi
    inizio = time.time()
    x, k = Jacobi(A,b,x0,toll,Kmax)
    fine = time.time()
    TempoJacobi[i] = fine - inizio
    IterazioniJacobi[i] = k
    
    # Calcolo tempo e iterazioni Gauss-Seidel
    inizio = time.time()
    x, k = GaussSeidel(A,b,x0,toll,Kmax)
    fine = time.time()
    TempoGaussSeidel[i] = fine - inizio
    IterazioniGaussSeidel[i] = k
    
    # Calcolo tempo MEG
    inizio = time.time()
    x = SoluzioneConMEG(A,b)
    fine = time.time()
    TempoMEG[i] = fine - inizio
    
    i = i + 1
    
# Rappresentazione dei tempi di calcolo
plt.figure(1)
plt.semilogy(nrange,TempoJacobi,label='Jacobi')
plt.semilogy(nrange,TempoGaussSeidel,label='Gauss-Seidel')
plt.semilogy(nrange,TempoMEG,label='MEG')
plt.xlabel('n')
plt.ylabel('Tempo di calcolo')
plt.legend()
plt.title('Confronto prestazioni')

# Rappresentazione del numero di iterazioni
plt.figure(2)
plt.plot(nrange,IterazioniJacobi,label='Jacobi')
plt.plot(nrange,IterazioniGaussSeidel,label='Gauss-Seidel')
plt.xlabel('n')
plt.ylabel('Iterazioni')
plt.legend()
plt.title('Confronto del numero di iterazioni')