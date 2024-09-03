import numpy as np
import time
import matplotlib.pyplot as plt
from sostituzione_avanti import SostituzioneAvantiOttimizzato
from sostituzione_indietro import SostituzioneIndietroOttimizzato

#Fattorizzazione LU
#==================
#Parametri di input:
#A :    matrice dei coefficienti
#==================
#Parametri di output:
#L :    parte triangolare inferiore
#U :    parte triangolare superiore

def Fatt_LU(A):
    A = np.copy(A)
    n = len(A)
    if abs(np.prod(np.diag(A))) < 1.0e-14:
        print('Errore: matrice possibilmente singolare!')
    else:
        for j in range(n-1):
            for i in range(j+1,n):
                A[i,j] = A[i,j]/A[j,j]
                for k in range(j+1,n):
                    A[i,k] = A[i,k] - A[i,j] * A[j,k]         
    L = np.tril(A,-1) + np.eye(n,n)
    U = np.triu(A)
    return L, U

def Fatt_LU_pivot(A):
	A = np.copy(A)
	n = len(A)
	indice = np.array(range(n))
	for j in range(n-1):
		A_max = abs(A[j,j]) ; i_pivot = j
		for i in range(j+1,n):
			if abs(A[i,j]) > A_max:
				A_max = abs(A[i,j]) ; i_pivot = i
		if i_pivot > j:
			for k in range(n):
				temp = A[i_pivot,k]
				A[i_pivot,k] = A[j,k]
				A[j,k] = temp
			temp = indice[i_pivot]
			indice[i_pivot] = indice[j]
			indice[j] = temp
		for i in range(j+1,n):
			A[i,j] = A[i,j]/A[j,j]
			for k in range(j+1,n):
				A[i,k] = A[i,k] - A[i,j]*A[j,k]
	L = np.tril(A,-1) + np.eye(n,n)
	U = np.triu(A)
	return L, U, indice

# Confronto delle prestazioni con matrice ad elementi casuali

def SoluzioneConLU(A,b):
    L, U = Fatt_LU(A)
    y = SostituzioneAvantiOttimizzato(L,b)
    x = SostituzioneIndietroOttimizzato(U,y)
    return x

def SoluzioneConLUPivot(A,b):
    L, U, i = Fatt_LU_pivot(A)
    y = SostituzioneAvantiOttimizzato(L,b[i])
    x = SostituzioneIndietroOttimizzato(U,y)
    return x

nmax = 100
nrange = range(5,nmax,10)
TempoLU = np.zeros(len(nrange))
TempoLUPivot = np.zeros(len(nrange))
ErroreLU = np.zeros(len(nrange))
ErroreLUPivot = np.zeros(len(nrange))
i = 0
for n in nrange:
    # Costruzione di un problema test
    A = (2*np.random.random((n,n))-1)*10
    xsol = np.ones((n,1))
    b = np.dot(A,xsol)
    
    # Calcolo tempo ed errore per fattorizzazione LU senza pivot
    itmax = 5
    itrange = range(itmax)
    Tempi = np.zeros(len(itrange))
    k = 0
    for it in itrange:  
        inizio = time.time()
        x = SoluzioneConLU(A,b)
        fine = time.time()
        Tempi[k] = fine - inizio
        k = k + 1
    TempoLU[i] = np.average(Tempi)
    ErroreLU[i] = np.linalg.norm(xsol-x)/np.linalg.norm(xsol)

    # Calcolo tempo ed errore per fattorizzazione LU con pivot
    Tempi = np.zeros(len(itrange))
    k = 0
    for it in itrange:  
        inizio = time.time()
        x = SoluzioneConLUPivot(A,b)
        fine = time.time()
        Tempi[k] = fine - inizio
        k = k + 1
    TempoLUPivot[i] = np.average(Tempi)
    ErroreLUPivot[i] = np.linalg.norm(xsol-x)/np.linalg.norm(xsol)

    i = i + 1
    
# Rappresentazione grafica del tempo di calcolo
plt.figure(1)
plt.semilogy(nrange,TempoLU,'r-',label='Senza pivot')
plt.semilogy(nrange,TempoLUPivot,'g-',label='Con pivot')
plt.xlabel('Dimensione del problema test')
plt.ylabel('Tempo di calcolo')
plt.legend()
plt.title('Confronto prestazioni fattorizzazione A=LU')
plt.show()

# Rappresentazione grafica dell'errore
plt.figure(2)
plt.semilogy(nrange,ErroreLU,'r-',label='Senza pivot')
plt.semilogy(nrange,ErroreLUPivot,'g-',label='Con pivot')
plt.xlabel('Dimensione del problema test')
plt.ylabel('Errore')
plt.legend()
plt.title('Errore al variare delle dimensioni')
plt.show()

# Confronto con matrice di Hilbert

nmax = 10
nrange = range(3,nmax)
TempoLU = np.zeros(len(nrange))
TempoLUPivot = np.zeros(len(nrange))
ErroreLU = np.zeros(len(nrange))
ErroreLUPivot = np.zeros(len(nrange))
i = 0
for n in nrange:
    # Costruzione di un problema test
    H = scipy.linalg.hilbert(n)
    xsol = np.ones((n,1))
    c = np.dot(H,xsol)
    
    # Calcolo tempo ed errore per fattorizzazione LU senza pivot
    itmax = 5
    itrange = range(itmax)
    Tempi = np.zeros(len(itrange))
    k = 0
    for it in itrange:  
        inizio = time.time()
        x = SoluzioneConLU(H,c)
        fine = time.time()
        Tempi[k] = fine - inizio
        k = k + 1
    TempoLU[i] = np.average(Tempi)
    ErroreLU[i] = np.linalg.norm(xsol-x)/np.linalg.norm(xsol)

    # Calcolo tempo ed errore per fattorizzazione LU con pivot
    Tempi = np.zeros(len(itrange))
    k = 0
    for it in itrange:  
        inizio = time.time()
        x = SoluzioneConLUPivot(H,c)
        fine = time.time()
        Tempi[k] = fine - inizio
        k = k + 1
    TempoLUPivot[i] = np.average(Tempi)
    ErroreLUPivot[i] = np.linalg.norm(xsol-x)/np.linalg.norm(xsol)

    i = i + 1
    
# Rappresentazione grafica del tempo di calcolo
plt.figure(1)
plt.semilogy(nrange,TempoLU,'r-',label='Senza pivot')
plt.semilogy(nrange,TempoLUPivot,'g-',label='Con pivot')
plt.xlabel('Dimensione del problema test')
plt.ylabel('Tempo di calcolo')
plt.legend()
plt.title('Confronto prestazioni fattorizzazione A=LU')
plt.show()

# Rappresentazione grafica dell'errore
plt.figure(2)
plt.semilogy(nrange,ErroreLU,'r-',label='Senza pivot')
plt.semilogy(nrange,ErroreLUPivot,'g-',label='Con pivot')
plt.xlabel('Dimensione del problema test')
plt.ylabel('Errore')
plt.legend()
plt.title('Errore al variare delle dimensioni')
plt.show()