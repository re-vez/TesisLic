#José Fernando Méndez Torres
#Códigos correspondientes al quinto capítulo de la tesis  
#"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.

#Este código regresa una molécula de ADN a partir de una 
#lista de lecturas y el ordenamiento de éstas. Además compara 
#a la molécula secuenciadacon la original.

import numpy as np
from Cap3 import noErrorFragmentAssembly, levinshtein
#Para ejecutar el código, deje el archivo Cap3.py en la misma 
#carpeta que éste.

with open("OrderDNA.npy",'rb') as f:
#Obtenemos el tour generado con MatrixToTour.py 
#u OrderGen.py    
    t = np.load(f)

t = t.astype(int)

#En esta lista almacenaremos las lecturas
Fragments = []

with open("Treated.txt","r") as f:
#Recordemos que por cada línea en el archivo hay una lectura
    for line in f:
        line = line.replace('\n', '')
        Fragments.append(line)

#Obtenemos la molécula a partir del orden de 
S = noErrorFragmentAssembly(Fragments,t)

#Extaemos la molécula original para compararlas
with open("Escherichia_phage_Lambda.txt","r") as f:
    _ = f.readline()
    DNA = ''
    for line in f:
        line = line.replace('\n', '')
        DNA += line

#Recordemos que necesitamos las longitudes para la proporción
lS = len(S)	
lDNA = len(DNA)

#Calculamos la distancia de Levenshtein
L = levinshtein(DNA,S)
#Calculamos la proporción
p = (1-L/max(lS,lDNA))
print("Hay un parecido de "+str(p*100)+" %.")

#Guardamos la molécula secuenciada
with open("Escherichia_phage_Lambda_seq_4.txt","w") as f:
    f.write("La distancia de Levenshtein es de "+str(L)+" pdb.\n")
    f.write("Hay un parecido del "+str(p*100)+"%.\n")
    f.write(S)


