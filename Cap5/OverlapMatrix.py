#José Fernando Méndez Torres
#Códigos correspondientes al quinto capítulo de la tesis  
#"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.
import numpy as np
from Cap3 import overlapMatrix
#Para ejecutar el código, deje el archivo Cap3.py en la misma 
#carpeta que éste.

#Aquí almacenaremos todas las lecturas del archivo .fastq
Fragments = []


with open("Treated.txt","r") as f:
    for line in f:
        line = line.replace('\n', '')
        Fragments.append(line)

C = overlapMatrix(Fragments)

with open("test.npy",'wb') as f:
    np.save(f, C)


