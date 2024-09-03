import numpy as np
import matplotlib.pyplot as plt

def TrapComposto(f,a,b,N):
    x = np.linspace(a,b,N+1)
    fx = f(x)
    S = 0
    for i in range(1,N):
        S = S + fx[i]
    TN = (b-a)/2/N*(fx[0] + 2*S + fx[N])
    return TN

def CavalieriSimpson(f,a,b,N):
    x = np.linspace(a,b,N+1)
    fx = f(x)
    S = 0
    for i in range(1,N):
        S = S + fx[i]
    SS = 0
    for i in range(1,N):
        SS = SS + f((x[i]+x[i+1])/2)
    SN = (b-a)/6/N*(f(a) + 2*S + 4*SS + f(b))
    return SN

# Funzione test da integrare
def f(x):
    y = x**3 + 6*x**2 - 7
    return y

# Primitiva della funzione test
def F(x):
    y = x**4/4 + 2*x**3 -7*x
    return y

# Intervallo di integrazione
a = -5.0 ; b = -2.0 ; 

# Calcolo integrale esatto
I = F(b) - F(a)

# Applicazione formule al variare di N
Nmax = 400
N_range = range(5,Nmax,10)
Err_Trap = np.zeros(len(N_range))
Err_Simp = np.zeros(len(N_range))
k = 0 
for N in N_range:
    TN = TrapComposto(f,a,b,N)
    SN = CavalieriSimpson(f,a,b,N)
    Err_Trap[k] = abs(I-TN)
    Err_Simp[k] = abs(I-SN)
    k = k + 1

# Rappresentazione grafica dell'errore al variare di N
plt.figure(1)
plt.semilogy(N_range,Err_Trap,'k-*',label='Trap. Comp.')
plt.semilogy(N_range,Err_Simp,'g-o',label='Cavalieri-Simpson')
plt.xlabel('N')
plt.ylabel('Errore')
plt.legend()
plt.title('Errore al variare di N')


