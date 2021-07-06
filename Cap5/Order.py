#José Fernando Méndez Torres
#Códigos correspondientes al quinto capítulo de la tesis  
#"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.

#Este código recibe recibe instancias del PAV codificadas en 
#una matriz de costos. Dichas instancias son resueltas con el 
#algoritmo de ramificación y acotamiento descrito en el 
#capítulo 2

import numpy as np
from BranchNBound import BnB

#Obtenemos la matriz guardada en un archivo .npy
with open("testTSP.npy",'rb') as f:
    C = np.load(f)

#Resolvemos la instancia
X,z = BnB(C)

#Guardamos el tour codificado en una matriz de adyacencia 
with open("order.npy",'wb') as f:
    np.save(f, X)
    