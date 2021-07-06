#José Fernando Méndez Torres
#Códigos correspondientes al quinto capítulo de la tesis  
#"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.

#Este código recibe recibe instancias del PAV codificadas en 
#una matriz de costos. Dichas instancias son resueltas con el 
#algoritmo de genético descrito en el capítulo 4

import numpy as np
from Genetics import geneticAlgorithmPlot
#Para ejecutar el código, deje el archivo BranchNBound.py en #la misma carpeta que éste.

#Obtenemos la matriz guardada en un archivo .npy
with open("test.npy",'rb') as f:
    C = np.load(f)
n = len(C)

#Encontramos un tour para la instancia
t, v = geneticAlgorithmPlot(C, 200, 10, 0.0005, 5000)

#Guardamos el tour codificado en una matriz de adyacencia 
with open("orderGen.npy",'wb') as f:
    np.save(f, t)


#El resto del código fue usado para crear la figura 5.5 
stringAux = ""

for i in range(len(v)):
    stringAux += "("+str(i)+","+str(v[i])+")"

print(stringAux)
with open("pgfplots.txt",'w') as f:
    f.write(stringAux)