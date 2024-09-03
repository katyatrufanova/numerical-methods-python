import numpy as np
import random
import matplotlib.pyplot as plt

# Dati di input
n = 200 # Dimensione dei vettori di input
x = np.array(random.sample(range(1, 10000),n))
x.sort()
y = np.array(random.sample(range(1, 10000),n))
y.sort()

# Costruzione dei dati perturbati
from random import random
xt = x + x*random()*1.0e-6 
yt = y + y*random()*1.0e-6

# Calcolo dei prodotti
P = np.multiply(x,y)
Pt = np.multiply(xt,yt)

# Calcolo errori sui dati e sui risultati
err_x = abs(x-xt)
err_y = abs(y-yt)
err_P = abs(P-Pt)

#Rappresentazione grafica
plt.figure(1)
plt.semilogy(err_x+err_y, err_P, 'g-')
plt.xlabel('Somma dell\'errore su x e dell\'errore su y')
plt.ylabel('Errore sul prodotto P')
plt.title('Rapporto tra errore sui dati ed errore sul risultato')
plt.show()