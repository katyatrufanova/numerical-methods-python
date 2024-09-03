from math import ceil
from math import log2
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

def BisezioniSuccessive(fun,a,b,tol,xs):
    # Verifica presenza di radici in [a,b]
    fa = fun(a) ; fb = fun(b)
    if fa*fb > 0:
        print('Errore: radice in [%f,%f] non garantita' % (a,b))
    else:
        n = ceil(log2((b-a)/tol)) - 1
        Err = np.zeros(n+1)
        for k in range(n+1):
            c = (a+b)/2 ; fc = fun(c)
            Err[k] = abs(xs-c)
            if fa*fc < 0 :
                b = c ; fb = fc
            else:
                a = c ; fa = fc 
    return c, Err

def Newton(fun,dfun,x0,tol,itmax):
    fx0 = fun(x0)
    dfx0 = dfun(x0)
    it = 0
    stop = 0
    n = ceil(log2((b-a)/tol)) - 1
    Err = np.zeros(n+1)
    while (not(stop) and it < itmax):
        x1 = x0 - fx0 / dfx0
        Err[it] = abs(xs-x1)
        fx1 = fun(x1)
        stop = abs(fx1) + abs(x1-x0)/abs(x1) < tol/5
        it = it + 1
        if (not(stop)):
            x0 = x1
            fx0 = fx1
            dfx0 = dfun(x0)
    if (not(stop)):
        print('Accuratezza richiesta non raggiunta in %d iterazioni ' % it)
    return x1, Err

def Secanti(fun,x1,x2,tol,itmax):
    fx1 = fun(x1)
    fx2 = fun(x2)
    it = 0
    stop = 0
    n = ceil(log2((b-a)/tol)) - 1
    Err = np.zeros(n+1)
    while (not(stop) and it < itmax):
        x3 = x2 - fx2 / ((fx2-fx1)/(x2-x1))
        Err[it] = abs(xs-x1)
        fx3 = fun(x3)
        stop = abs(fx3) + abs(x3-x2)/abs(x3) < tol/5
        it = it + 1
        if (not(stop)):
            x1 = x2
            x2 = x3
            fx1 = fx2
            fx2 = fx3
    if (not(stop)):
        print('Accuratezza richiesta non raggiunta in %d iterazioni ' % it)
    return x3, Err

def Corde(fun,x0,m,tol,itmax):
    fx0 = fun(x0)
    it = 0
    stop = 0
    n = ceil(log2((b-a)/tol)) - 1
    Err = np.zeros(n+1)
    while (not(stop) and it < itmax):
        x1 = x0 - fx0 / m
        Err[it] = abs(xs-x1)
        fx1 = fun(x1)
        stop = abs(fx1) + abs(x1-x0)/abs(x1) < tol/5
        it = it + 1
        if (not(stop)):
            x0 = x1
            fx0 = fx1
    if (not(stop)):
        print('Accuratezza richiesta non raggiunta in %d iterazioni ' % it)
    return x1, Err

def fun(x):
    y = x**2 - 2
    return y

def dfun(x):
    y = 2*x
    return y

xs = sqrt(2)
x0 = float(input('Inserire l\'approssimazione iniziale x0: ')) 
x2 = float(input('Inserire l\'approssimazione iniziale x2 per il metodo delle Secanti: ')) 
m = float(input('Inserire il coefficiente angolare m per il metodo delle corde: ')) 

itmax = 1000
a = 0 ; b = 2 ; tol = 1.0e-8
x_BS, Err_BS = BisezioniSuccessive(fun,a,b,tol,xs)
x_N, Err_N = Newton(fun,dfun,x0,tol,itmax)
x_S, Err_S = Secanti(fun,x0,x2,tol,itmax)
x_C, Err_C = Corde(fun,x0,m,tol,itmax)

# Rapprentazione grafica dell'errore
plt.figure(1)
plt.semilogy(range(len(Err_BS)),Err_BS,label='Errore Bisezioni Successive')
plt.semilogy(range(len(Err_N)),Err_N,label='Errore Newton')
plt.semilogy(range(len(Err_S)),Err_S,label='Errore Secanti')
plt.semilogy(range(len(Err_C)),Err_C,label='Errore Corde')
plt.semilogy(np.array([0,len(Err_BS)]),np.array([tol,tol]),'r:')
plt.title('Confronto degli errori')
plt.xlabel('n')
plt.legend()
plt.show()