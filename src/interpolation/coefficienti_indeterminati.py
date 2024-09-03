import numpy as np
import matplotlib.pyplot as plt
from metodi_iterativi import SoluzioneConMEG

def fun(x):
    y = 2*x**5 + np.sin(x) + np.cos(x)    
    return y

# Intervallo di interpolazione
a = -2 ; b = 2

# Grado ed i nodi di interpolazione
n = 20
xn = np.linspace(a,b,n+1)
yn = fun(xn)

# Valori in cui calcolare la funzione e il polinomio di interpolazione
x = np.linspace(a,b,200)
fx = fun(x)

# Calcolo polinomio mediante metodo dei coefficienti indeterminati
Vn = np.vander(xn)
c = SoluzioneConMEG(Vn,yn)
px = np.polyval(c,x)

# Grafico del polinomio di interpolazione e della funzione
plt.figure(1)
plt.plot(x,fx,'b-',label='f(x)')
plt.plot(x,px,'r--',label='p_n(x)')
plt.plot(xn,yn,'ro')
plt.xlabel('x')
plt.legend()
plt.show()


