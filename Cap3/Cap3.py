#José Fernando Méndez Torres
#Códigos correspondientes al tercer capítulo de la tesis  
#"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.

import numpy as np

def levinshtein(u, v):
    """
    La función levenshtein recibe dos cadenas de caracteres
    u,v y regresa la distancia de levenshtein d(u,v). Ésta es 
    definida como la mínima cantidad de cambios que se debe
    realizar a u para obtener a v, donde estos cambios pueden 
    ser la inclusión, sustitución o eliminación de un 
    caracter. Para realizar el cálculo, se utiliza 
    programación dinámica. 
    """
    n = len(u)+1
    m = len(v)+1
    d = np.zeros((n,m))
    #La matriz d nos ayudará a almacenar cálculos, pues la 
    #distancia se computa de manera recursiva
    for i in range(n):
        d[i][0]=i
    for i in range(m):
        d[0][i] = i
    #Los ciclos for's anteriores son la base recursiva del 
    #cálculo
    for i in range(1, n):
        for j in range(1, m):
            #Ecuación recursiva
            d[i][j] = min(d[i][j-1]+1, d[i-1][j]+1, d[i-1][j-1]+(not u[i-1] == v[j-1]))
    #La última entrada de la matriz es igual a d(u,v)
    return float(d[n-1][m-1])

def overlap(a, b, length=0):
    """
    Dadas dos cadenas u,v, esta función nos regresa la 
    longitud del traslape máximo entre ambas cadenas. 
    Adicionalmente podemos solicitar una longitud mínima de 
    traslape y en caso de no encontrar el traslape, el método 
    nos regresa como resultado 0.
    """
    aux = 0
    while True:
    #La función permanecerá en el ciclo hasta satisfacer 
    #alguno de los if's
        aux = a.find(b[:length],aux)
        #find nos indica el índice en el que se encuentra la 
        #primera coincidencia con el prefijo de b de longitud 
        #length a partir del índice aux 
        if aux == -1:
        #En caso de que no existan coincidencias, el método 
        #find arrojará el número -1, esto indica que el 
        #traslape es la cadena vacía
            return 0
        if b.startswith(a[aux:]):
        #Si el sufijo que empieza en el índice aux de la 
        #cadena a resulta ser prefijo de b, regresamos |a|-aux
            return len(a)-aux
        aux += 1



def overlapMatrix(reads, k = 0):
    """
    Dada una lista de lecturas y una longitud mínima de 
    traslapes, esta función transforma el problema de 
    secuenciación de novo en una instancia del problema del 
    mensajero. Al resolver dicho problema, podemos obtener 
    una secuencia candidata correspondiente a las lecturas. 
    """
    n = len(reads)
    C = np.array(n*[n*[np.infty]])
    for i in range(n):
        for j in range(n):
            if i != j:
            #La diagonal va con costos infinitos
                C[i,j] = len(reads[i])+len(reads[j])-2*overlap(reads[i],reads[j],k)
                #Recordemos que los costos deben ser no 
                #negativos y queremos minimizar la función 
                #objetivo, por lo que se define a los costos 
                #de la forma anterior
    return(C)
        
def noErrorFragmentAssembly(reads, indexes):
    """
    Dada una lista de lecturas y un ordenamiento de ellas 
    (preferentemente obtenido al resolver el problema del 
    mensajero con la matriz de traslapes), obtenemos la 
    secuenciación de las lecturas.
    """
    seq = reads[indexes[0]]
    #Empezamos con la primer cadena del ordenamiento
    for i in range(len(reads)-1):
        O = overlap(reads[indexes[i]],reads[indexes[i+1]],1)
        seq += reads[indexes[i+1]][O:]
        #Posteriormente, iremos agregando el resto 
        #concatenando con los sufijos del resto de las cadenas
    return seq

def acurratePercentage(u, v, L):
    """
    Dadas dos cadenas u,v y L la distancia de Levenshtein 
    entre ellas, podemos ver en qué proporción son semejantes 
    con p = 1-L/max{|u|,|v|}.
    """
    return 1-L/max(len(u),len(v))

def errorOverlapCost(u, v):
    """
    Dadas dos cadenas u,v, se puede calcular una longitud 
    penalizada del mejor traslape aproximado bajo el criterio:
    max{|p|+|s|-k*d(p,s): p es sufijo de u, s es prefijo de v}
    con k = (|u|+|v|)/max{|u|,|v|}. Esta rutina puede ser 
    utilizada para un algoritmo de secuenciación de novo. 
    Además, la función regresa a los mejores candidatos a 
    sufijo de u y prefijo de v.
    """
    cost = 0
    #Al querer maximizar, empezamos con un costo cero.
    n = len(u)
    m = len(v)
    s = '' 
    p = ''
    #Nos apoyaremos de las cadenas vacías anteriores para 
    #almacenar a los mejores candidatos a traslapes 
    #aproximados
    for i in range(n):
        for j in range(1,m+1):
        #utilizando estos ciclos for's, buscaremos sobre 
        #todas las subcadenas posibles
            su = u[i:n+1]
            pv = v[0:j]
            k = (j+n-i)/max(j,n-i)
            c = j+(n-i)-k*levinshtein(su,pv)
            #Al mejorar el costo, actualizamos cost y las 
            #cadenas s,p
            if c > cost:
                cost = c
                s = su
                p = pv
    print(cost)
    finalCost = m+n-cost 
    #Realizamos este cambio pues queremos resolver un 
    #problema de minimización (PAV)
    return finalCost,s,p
