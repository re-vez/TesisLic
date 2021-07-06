#José Fernando Méndez Torres
#Códigos correspondientes al quinto capítulo de la tesis  
#"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.

#Este archivo utiliza una lista de lecturas tratadas 
#(revise NoErrorTreatment.py) para generar una matriz que 
#codifique una instancia del problema del mensajero

import numpy as np
from Cap3 import overlapMatrix
#Para ejecutar el código, deje el archivo Cap3.py en la misma 
#carpeta que éste.

#Aquí almacenaremos todas las lecturas del archivo .txt
Fragments = []

with open("Treated.txt","r") as f:
    for line in f:
        line = line.replace('\n', '')
        Fragments.append(line)
        #Las vamos recolectando línea por línea

C = overlapMatrix(Fragments)
#Obtenemos la matriz de traslapes

with open("test.npy",'wb') as f:
#La guardamos en el archivo test.py
    np.save(f, C)


