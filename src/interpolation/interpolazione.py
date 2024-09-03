import numpy as np
import matplotlib.pyplot as plt

# Interpolazione di Lagrange mediante formula baricentrica
# ========================================================
# Parametri di input
# xn : vettore dei nodi
# yn : vettore dei valori associati ai nodi
# x : vettore di punti in cui calcolare il polinomi di interpolazione
# ========================================================
# Parametri di output
# p : vettore di punti con il valore del polinomio in x
# ========================================================
def Lagrange(xn,yn,x):
    
    # Calcolo coefficienti formula baricentrica
    zn = Z_Coeff(xn,yn)
    
    # Calcolo polinomio interpolazione nei punti di x
    p = np.zeros(len(x))
    for i in range(len(x)):
        p[i] = CalcLagrange(x[i],xn,yn,zn)
        
    return p

# Calcolo dei coefficienti formula baricentria
def Z_Coeff(xn,yn):
    n = len(xn)
    X = np.eye(n)
    for i in range(n):
        for j in range(n):
            if j > i:
                X[i,j] = xn[i]-xn[j]
            elif j < i:
                X[i,j] = - X[j,i]
    zn = np.zeros(n)
    for j in range(n):
        zn[j] = yn[j]/np.prod(X[j,:])
    return zn

# Calcolare il polinomio di interpolazione in un punto x
def CalcLagrange(x,xn,yn,zn):
    check_nodi = abs(x-xn) < 1.0e-14
    if True in check_nodi: 
        temp = np.flatnonzero(check_nodi == True)
        j = temp[0]
        pn = yn[j]
    else:
        n = len(xn)
        S = 0
        for j in range(n):
            S = S + zn[j]/(x-xn[j])
        pn = np.prod(x-xn)*S
    return pn

# Calcolo delle differenze divise,
# ossia i coefficienti del polinomio di Newton  
# ===============================
# Parametri di input:
# ===============================
# xn : vettore dei nodi di interpolazione
# yn : vettore dei valori associati
# ===============================
# Parametri di output:
# ===============================
# vettore dei coefficienti del polinomio di Newton (differenze divise) 

def CoefficientiNewton(xn,yn):
    n = len(xn)
    A = np.zeros((n,n))
    for i in range(0,n):
        A[i,0] = yn[i]
    for j in range(1,n+1):
         for i in range(j,n):
             A[i,j] = (A[i,j-1] - A[i-1,j-1]) / (xn[i] - xn[i-j])
    return np.diagonal(A)

# Calcolo del polinomio di Newton
# ===============================
# Parametri di input:
# ===============================
# x : punto in cui calcolare il polinomio
# xn : vettore dei nodi di interpolazione
# d : vettore dei coefficienti
# ===============================
# Parametri di output:
# ===============================
# p : polinomio di Newton calcolato nel punto x

def CalcoloPolinomioNewton(x, xn, d):
    n = len(d)
    p = d[-1]
    for i in range(n-1,-1,-1):
        p = d[i] + p*(x-xn[i])
    return p

def Newton(xn,yn,x):
    d = CoefficientiNewton(xn,yn)
    px = CalcoloPolinomioNewton(x,xn,d)
    return px

def funz(x):
    y = x**2 + np.cos(2*x)
    return y

# Dati di interpolazione
a = -np.pi ; b = np.pi
n = 30 # Grado del polinomio
xn = np.linspace(a,b,n+1)
yn = funz(xn)

# Calcolo funzione e polinomio di interpolazione
x = np.linspace(a,b,200)
px = Newton(xn,yn,x)
fx = funz(x)

# Grafico funzione e polinomio
plt.figure(1)
plt.plot(x,fx,'g-',label='f(x)')
plt.plot(x,px,'r--',label='pn(x)')
plt.plot(xn,yn,'ro')
plt.xlabel('x')
plt.legend()
plt.title('Polinomio di interpolazione di Newton')
plt.show()

# Grafico del resto
plt.figure(2)
plt.semilogy(x,abs(fx-px),label='|f(x)-pn(x)|')
plt.xlabel('x')
plt.legend()
plt.title('Resto di interpolazione')
plt.show()

# Confronto delle formule di interpolazione

def funz(x):
    y = x**2 + np.cos(2*x)
    return y

# Dati di interpolazione
a = -np.pi ; b = np.pi

# Calcolo della funzione
x = np.linspace(a,b,400)
fx = funz(x)

# Calcolo dei polinomi e resto al variare di n con nodi equidistanti
nmax = 60
RL = np.zeros(nmax)
RN = np.zeros(nmax)
for n in range(nmax):
    
    # Calcolo nodi e valore funzione nei nodi
    xn = np.linspace(a,b,n+1)
    yn = funz(xn)
    
    # Calcolo polinomi di interpolazione di grado n
    pxL = Lagrange(xn,yn,x)
    pxN = Newton(xn,yn,x)
    
    # Calcolo max Resto per pol. interp. grado n
    RL[n] = max(abs(fx-pxL))
    RN[n] = max(abs(fx-pxN))
    
# Rappresentazione grafica del resto al variare di n con nodi equidistanti
plt.figure(1)
plt.semilogy(range(nmax),RL,label='Resto con Lagrange')
plt.semilogy(range(nmax),RN,label='Resto con Newton')
plt.xlabel('n')
plt.legend()
plt.title('Resto con nodi equidistanti')
plt.show()

# Calcolo dei polinomi e resto al variare di n con nodi di Chebyshev
nmax = 60
RL = np.zeros(nmax)
RN = np.zeros(nmax)
for n in range(nmax):
    
    # Calcolo nodi e valore funzione nei nodi
    k = np.array(range(n,-1,-1))
    t = np.cos((2*k+1)*np.pi/2/(n+1)) 
    xn = (a+b)/2 + (b-a)/2*t
    yn = funz(xn)
    
    # Calcolo polinomi di interpolazione di grado n
    pxL = Lagrange(xn,yn,x)
    pxN = Newton(xn,yn,x)
    
    # Calcolo max Resto per pol. interp. grado n
    RL[n] = max(abs(fx-pxL))
    RN[n] = max(abs(fx-pxN))
    
# Rappresentazione grafica del resto al variare di n con nodi di Chebyshev
plt.figure(2)
plt.semilogy(range(nmax),RL,label='Resto con Lagrange')
plt.semilogy(range(nmax),RN,label='Resto con Newton')
plt.xlabel('n')
plt.legend()
plt.title('Resto con nodi di Chebyshev')
plt.show()

