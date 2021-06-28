#José Fernando Méndez Torres
#Códigos correspondientes al cuarto capítulo de la tesis  
#"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.

import numpy as np 
import matplotlib.pyplot as plt

def random_chromosome(nCities):
    """
    Esta función recibe un número de ciudades y regresa un 
    camino hamiltoniano aleatorio de ellas.
    """
    path = np.arange(nCities) #El recorrido 1,2,...,n
    np.random.shuffle(path)     #Ordenamiento aleatorio
    return np.array(path)

def initial_population(popSize,nCities):
    """
    Esta función recibe una cantidad de caminos popSize y un 
    número de ciudades para regresar popSize caminos 
    aleatorios para recorrer las n ciudades
    """
    population = []
    for _ in range(popSize):
        population.append(random_chromosome(nCities))
    return np.array(population)

def path_eval(Path, CostMatrix):
    """
    Esta función nos regresa el costo de un camino dada una
    matriz de costos de ir de una ciudad i a una ciudad j.
    """
    nCities = len(Path) #El número de ciudades a visitar
    cost = 0 #En esta variable almacenaremos el costo
    for i in range(nCities-1):
        cost += CostMatrix[Path[i],Path[i+1]]
    return float(cost)

def path_fitness(tour,costMatrix):
    """
    Dado un camino y una matriz de costos, esta función nos 
    regresa la aptitud del camino. Recordemos que a mayor 
    aptitud, mayor es la probabilidad de elegir a un poblador 
    para su cruza y que el PM es un problema de minimización, 
    por ello regresamos al inverso del costo del camino
    """
    fitness = 1/path_eval(tour,costMatrix)
    return fitness

def pop_fit(Pop, CostMatrix):
    """
    Dada una tupla con caminos Pop y una matriz de costos, 
    esta función nos regresa un vector de aptitudes v y 
    además ordenamos a v y Pop con respecto a v.
    """
    nPop = len(Pop) #Obtenemos la cantidad de pobladores
    Apt = np.zeros(nPop) #Inicializamos el vector de aptitudes
    for i in range(nPop):
    #Para cada poblador, obtenemos su aptitud
        Apt[i] = path_fitness(Pop[i], CostMatrix)
    #Obtenemos los índices de reordenamiento a partir de Apt
    Indexes = np.argsort(Apt)[::-1]
    #Reordenamos a la población y a Apt
    return Pop[Indexes], Apt[Indexes]

def selection_index(Apt, elite):
    """
    Dado un vector de aptitudes y un parámetro de elitismo, 
    vamos a obtener m-elite índices para generar nuevos 
    cromosomas, donde m es la dimensión de Apt.
    """
    #Construimos una función de distribución a partir de Apt
    F = np.zeros(len(Apt)) 
    S = sum(Apt)
    F[0] = Apt[1]/S
    for i in range(1, len(Apt)-1):
        F[i] = F[i-1] + Apt[i]/S
    F[len(Apt)-1] = 1
    #A partir de aquí, inicializamos el vector de índices
    nIndexes = len(Apt)-elite
    Indexes = np.array(nIndexes*[0])
    #Aplicamos lo descrito por G. Pagès 
    #(Teorema fundamental de la simulación)
    for i in range(nIndexes):
        eta = 0
        xi = np.random.random()
        while xi > F[eta]:
            eta += 1
        #Obtenemos cada índice a partir de números 
        #pseudoaleatorios
        Indexes[i] = eta
    return Indexes

def breed_mpx(Pop1, Pop2):
    """
    A partir de dos cromosomas, podemos generar un nuevo 
    cromosoma con el cruce conservativo maximal (MPX). 
    """
    nCities = len(Pop1) #Obtenemos el número de ciudades
    NewPop = [] #Aquí almacenaremos el inicio del camino
    #Los siguientes números aleatorios servirán para tomar un 
    #tramo del primer camino, desde el índice m hasta el M-1. 
    xi1 = np.random.randint(0,nCities-1) 
    xi2 = np.random.randint(0,nCities-1)
    m = min(xi1,xi2)
    M = max(xi1,xi2)
    #Recuperamos la info de Pop1
    for i in range(m,M):
        NewPop.append(Pop1[i])
    #Recorremos las ciudades faltantes en el orden de Pop2
    N = [item for item in Pop2 if item not in NewPop]
    return np.array(NewPop + N)
    
def breed_pop(Pop, Ind, elite):
    """
    Esta función genera una nueva población a partir de una 
    población Pop, un conjunto de índices de selección Ind y 
    un parámetro de elitismo.
    """
    nPop = len(Pop) #Obtenemos el número de pobladores
    NewPop = []     #Aquí almacenaremos a cada poblador
    for i in range(elite):
    #Conservamos a los mejores elite pobladores de Pop
        NewPop.append(Pop[i])
    for i in range(nPop - elite):
    #Generamos nPop-elite nuevos pobladores
        sigma = breed_mpx(Pop[Ind[i-1]],Pop[Ind[i]])
        NewPop.append(sigma)
    #Convertimos a la lista en un arreglo de numpy
    return np.array(NewPop)

def mutation(sigma):
    """
    Esta función modifica a un ordenamiento con por medio de 
    la mutación por intercambio: Dados dos ciudades 
    aleatorias, se intercambian de lugar.
    """
    nCities = len(sigma)
    #Obtenemos las dos ciudades aleatorias
    xi1 = np.random.randint(0,nCities-1) 
    xi2 = np.random.randint(0,nCities-1) 
    #Realizamos el intercambio
    h = sigma[xi1]
    sigma[xi1] = sigma[xi2]
    sigma[xi2] = h
    return sigma

def pop_mutation(Pop, pm):
    """
    Esta función muta a cada poblador con probabilidad pm.
    La mutación es por intercambio (de ciudades).
    """
    for i in range(len(Pop)):
        #Para cada poblador, generamos un número aleatorio
        if np.random.random() <= pm:
        #Si este resulta ser menor o igual a pm, mutamos al 
        #i-ésimo poblador
            Pop[i] = mutation(Pop[i])
    return Pop

def next_generation(Pop, CostMatrix, elite, pm):
    """
    Esta función genera una nueva población con una población 
    Pop, una matriz de costos CostMatrix, un parámetro de 
    elitismo elite y una probabilidad de mutación pm.
    """
    #Primero generamos un vector de aptitudes y ordenamos a 
    #ambos vectores de manera descendente
    Pop, Apt = pop_fit(Pop, CostMatrix)
    #Obtenemos un vector de índices m-elite índices que 
    #utilizaremos para seleccionar a los padres de la próxima 
    #generación, con m la dimensión de Apt
    Indexes = selection_index(Apt, elite)
    #Generamos m-elite nuevos pobladores y conservamos a los 
    #elite mejores
    Pop = breed_pop(Pop, Indexes, elite)
    #Sometemos a la población a un proceso de mutación
    Pop = pop_mutation(Pop, pm)
    return Pop

def genetic_algorithm(PopSize, CostMatrix, elite, pm, gen):
    """
    Este algoritmo genético obtiene 'buenas' soluciones para 
    instancias del PM codificadas en una matriz de costos 
    CostMatrix. Además de dicha matriz, los otros parámetros 
    que recibe son el tamaño de la población PopSize, la 
    probabilidad de mutación pm y el número de generaciones 
    a simular. Como resultado da la mejor instancia 
    encontrada y su valor.
    """
    #Obtenemos el número de ciudades y generamos a la
    #población inicial
    nCities = len(CostMatrix[0])
    Pop = initial_population(PopSize, nCities)
    #En los siguientes auxiliares guardaremos la mejor 
    #solución encontrada y su mejor valor
    bestPath = np.array(nCities*[0])
    bestValue = np.inf
    #A continuación, simularemos la dinámica poblacional
    for i in range(0, gen):
        Pop = next_generation(Pop, CostMatrix, elite, pm)
        Pop, Apt = pop_fit(Pop, CostMatrix)
        #Evaluamos a la población y verificamos si el mejor
        #cromosoma supera a la mejor solución encontrada
        #hasta el momento
        if 1/Apt[0] < bestValue:
            bestPath = Pop[0]
            bestValue = 1/Apt[0]
        #Regresamos la mejor solución encontrada y su valor
    return bestPath, bestValue



def genetic_algorithm_plot(PopSize, CostMatrix, elite, pm, gen):
    """
    Este algoritmo genético obtiene 'buenas' soluciones para 
    instancias del PM codificadas en una matriz de costos 
    CostMatrix. Además de dicha matriz, los otros parámetros 
    que recibe son el tamaño de la población PopSize, la 
    probabilidad de mutación pm y el número de generaciones 
    a simular. Como resultado da la mejor instancia 
    encontrada y su valor. Además genera una gráfica de los 
    mejores valores de cada generación para visualizar el 
    desempeño del método.
    """
    #Obtenemos el número de ciudades y generamos a la
    #población inicial
    nCities = len(CostMatrix[0])
    Pop = initial_population(PopSize, nCities)
    #En los siguientes auxiliares guardaremos la mejor 
    #solución encontrada, su mejor valor y el progreso
    #de la dinámica poblacional
    bestPath = np.array(nCities*[0])
    bestValue = np.inf
    progress = []
    #A continuación, simularemos la dinámica poblacional
    for i in range(0, gen):
        Pop = next_generation(Pop, CostMatrix, elite, pm)
        Pop, Apt = pop_fit(Pop, CostMatrix)
        #Evaluamos a la población y verificamos si el mejor
        #cromosoma supera a la mejor solución encontrada
        #hasta el momento
        if 1/Apt[0] < bestValue:
            bestPath = Pop[0]
            bestValue = 1/Apt[0]
        #Almacenamos el valor del mejor cromosoma
        progress.append(1/Apt[0])
    #Finalmente generamos y ...
    plt.plot(progress)
    plt.ylabel('Distance')
    plt.xlabel('Generation')
    plt.show()
    #Regresamos la mejor solución encontrada y su valor
    return bestPath, bestValue

