#José Fernando Méndez Torres
#Códigos correspondientes al quinto capítulo de la tesis  
#"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.

#Este código regresa el camino de menor costo del PM asociado 
#al tour de menor costo del PM codificado como matriz
import numpy as np

#Obtenemos a la matriz resultante de correr Order.py
with open("orderTSP7.npy","rb") as f:
    rawX = np.load(f)

#La cantidad de lecturas
n = len(rawX)-1

#Fijamos al índice j = 0
j = 0
#Este contador nos ayudará a encontrar el extremo de la 
#caminata
count = np.inf

while count == np.inf:
    #Mientras no encontremos el extremo
	if rawX[n,j] == 1:
        #Para este punto, encontramos al extremo
		count = j
	else:
        #En otro caso, pasamos al sucesor de j
		j += 1


count = 0 #rederinimos al contador
t = np.zeros(n) #en esta variable almacenaremos la caminata
t[count] = j #Almacenamos al extremo de la caminata
i = j #definimos a i = j
j = 0 #y redefinimos a j
#Eliminamos la última fila y la última columna de la matriz, 
#pues ya no nos serán útiles
newX = np.delete(rawX, n, 0)
X = np.delete(newX, n, 1)

while count < n-1:
#Mientras no hayamos recorrido todas las ciudades
    while X[i, j] == 0:
    #si la conexión entre i y j no está en el tour
        j += 1
        #aumentamos en 1 al índice j
    else:
        #almacenamos la ciudad j
        t[count+1] = j
        #reiniciamos los índices
        i = j
        j = 0
        #aumentamos en 1 al contador de ciudades
        count += 1

#np.zeros() arroja un arreglo de flotantes, como t es un #arreglo de índices, necesitamos que sea de enteros
t = t.astype(int)

#Guardamos el tour 
with open("OrderDNA.npy",'wb') as f:
    np.save(f, t)
