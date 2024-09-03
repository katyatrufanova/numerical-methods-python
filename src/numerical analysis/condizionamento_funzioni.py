import numpy as np
import matplotlib.pyplot as plt

'''
Il seguente codice presenta lo studio del condizionamento delle funzioni x3, x(x2-1) e 1/(2-3x). 
Il condizionamento del calcolo di una funzione f(x) dipende sia dalla funzione che dai punti in cui viene calcolata, come mostrato dall’output del codice presentato.

L’indice di condizionamento della funzione x3 è uguale a |x * (f’(x) / f(x))| = |x * (3x2 / x3)| = 3,
dunque la funzione è sempre ben condizionata.

L’indice di condizionamento di x(x2-1) è uguale a |x * (f’(x) / f(x))| = |x * ( 3x2-1 / x(x2-1) )| =       
| 3x2-1 / x2-1 |, quindi la funzione è mal condizionata se x è vicino a 1 o a -1, mentre risulta ben condizionata per gli altri valori di x.

L’indice di condizionamento di 1/(2-3x) è uguale a |x * (f’(x) / f(x))| = |x * (3 / (2-3x)2) * (2-3x)| = |3x / 2-3|, 
quindi la funzione risulta mal condizionata quando x è vicino a 2/3 e ben condizionata altrimenti.
'''

def f1(x):
    return x**3

def df1(x):
    return 3*x**2

def f2(x):
    return x*(x**2-1)

def df2(x):
    return 3*x**2 - 1

def f3(x):
    return 1/(2-3*x)

def df3(x):
    return 3/(2-3*x)**2

def IndiceCondizionamento(x,f,df):
    return abs(x*(df(x)/f(x)))

x = np.linspace(-20,20,200)

# Calcolo degli indici di condizionamento
indCond_f1 = IndiceCondizionamento(x,f1,df1)
indCond_f2 = IndiceCondizionamento(x,f2,df2)
indCond_f3 = IndiceCondizionamento(x,f3,df3)

# Rappresentazione grafica dei risultati
plt.figure(1)
plt.plot(x, indCond_f1, label='$x^3$')
plt.plot(x, indCond_f2, label='x($x^2$-1)')
plt.plot(x, indCond_f3, label=r'$\frac{1}{2-3x}$')
plt.xlabel('x')
plt.ylabel('Indice di condizionamento')
plt.title('Condizionamento del calcolo del valore di una funzione')
plt.legend()
plt.show()