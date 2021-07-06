#José Fernando Méndez Torres
#Códigos correspondientes al quinto capítulo de la tesis  
#"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.

#Este programa transforma (en tiempo polinomial) matrices que 
#codifican una instancia del PM a una equivalente del PAV.

import numpy as np

#Obtenemos la instancia del PM codificada en una matriz
with open("test.npy","rb") as f:
    C = np.load(f)

n = len(C)
C1 = np.zeros((n, 1))
Cnew = np.hstack((C, C1))
#Aumentamos una fila (con ceros)
C2 = np.zeros(n+1)
Cnew = np.vstack([Cnew, C2])
#Aumentamos una columna (con ceros) e infinitos en la diagonal
Cnew[n, n] = np.inf


#Guardamos la instancia del PAV
with open("testTSP.npy",'wb') as f:
    np.save(f, Cnew)