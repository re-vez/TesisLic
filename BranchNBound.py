#José Fernando Méndez Torres
#Códigos correspondientes al segundo capítulo de la tesis  #"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.

import numpy as np
import copy
            
def Hungarian(C):
    """
    Este algoritmo resuelve problemas de asignación a partir 
    de un modelo lineal en tiempo O(n^3), donde n es la
    cantidad de tareas.
    Recordemos que un problema de asignación aquel en donde
    dadas n tareas y n personas para realizarlas, queremos
    repartir las tareas para efectuarlas de la manera más
    eficiente bajo las restricciones de que una sola persona
    puede resolver una tarea y todas ellas deben ser 
    efectuadas. Notemos que lo anterior puede ser representado
    como un acoplamiento perfecto en gráficas de menor costo.
    El algoritmo es de tipo primal-dual y consiste de tres
    subrutinas, en la primera se encuentra una solución 
    factible del problema dual y es llamada preprocesamiento,
    en la segunda se verifica si se han asignado todas las 
    tareas y es llamada húngaro y la tercera encuentra 
    cadenas o caminos aumentantes para introducir una tarea 
    no asignada y es llamada procedimiento alternante.
    """
    ##########################################################
    #Preprocesamiento
    n = len(C[0]) #La cantidad de tareas
    row = n*[np.inf]
    phi = n*[np.inf]

    #row y phi nos ayudarán a mantener el seguimiento de las
    #entradas X's, recordemos que row[j] = i iff phi[i] = j  
    #iff X[i,j]=1 iff se asigna la tarea i a la persona j

    u = n*[0]
    v = n*[0]
    #Las anteriores son las variables duales 

    V = set(range(n))
    Q = set([]) #El conjunto de tareas asignadas
    x = np.array(n*[n*[0]])
    for i in range(n):
        u[i] = min(C[i,:])
    for j in range(n):
        v[j] = min(C[:,j]-u[:])
    #Con lo anterior, u,v continúa siendo una solución dual
    #factible
    for i in range(n):
        for j in range(n):
            if row[j] == np.inf and C[i,j]-u[i]-v[j]==0 and(i not in Q):
            #En caso de que la j-ésima persona esté libre,
            #de que el costo reducido sea 0
            #y de que la tarea i aún no esté asignada,
            #asignamos la tarea i a la persona j.
                row[j] = i
                phi[i] = j
                Q = Q|set([i])

    ##########################################################
    #Húngaro 1/2
    while len(Q) < n:
        k = list(V-Q)[0] #dada una tarea sin asignar...
    ##########################################################
    #Procedimiento alternante
        pi = np.array(n*[np.inf])
        #El vector pi nos será de utilidad para actualizar
        #los vectores duales u,v
        SU = set([])  #tareas examinadas
        LV = set([])  #personas etiquetadas 
        SV = set([])  #personas examinadas 
        sink = np.inf 
        #sink es un auxiliar que utilizaremos para almacenar 
        #a una persona sin tarea asignada 
        i = k
        #k es la tarea elegida para asignar en este paso
        index = 1
        pred = n*[i]
        #El vector pred nos ayudará a actualizar el 
        #acoplamiento cuando encontremos una tarea sin asignar
        while sink == np.inf:
            SU = SU|set([i])
            #Marcamos a la tarea i como examinada
            for w in V - LV:
            #A todas las personas no examinadas
                if C[i,w]-u[i]-v[w]< pi[w]:
                # si el costo reducido de i,w es menor que
                #pi[w], entonces
                    pred[w] = i #redefinimos a pred[w]
                    pi[w] = C[i,w]-u[i]-v[w]
                    #actualizamos el costo reducido 
                    if pi[w] == 0:
                    #si pi[w] es cero, marcamos a la 
                    #persona como etiquetada
                        LV = LV|set([w])
            
            if LV - SV == set([]):
            #Si todas las personas etiquetadas han sido 
            #examinadas...
                delta = np.inf
                for w in V-LV:
                    if delta > pi[w]:
                        delta = pi[w]
                #Definimos a delta := min{pi[j]:j in V\LV}
                for w in SU:
                    u[w] = u[w]+delta
                for w in LV:
                    v[w] = v[w]-delta
                #Actualizamos los vectores duales
                for w in V-LV:
                    pi[w] = pi[w]-delta
                #Actualizamos a pi en las entradas no 
                #etiquetadas
                    if pi[w] == 0:
                    #Si el costo reducido es cero
                    #etiquetamos a la persona w
                        LV = LV|set([w])
            j = list(LV-SV)[0] 
            #Sea j una persona etiquetada y no examinada
            SV = SV|set([j])
            #La marcamos como examinada
            if row[j] == np.inf:
            #Si j no tiene tareas asignadas, asignamos a sink
                sink = j
            else:
                i = row[j]
            index += 1    
    ##########################################################
    # Húngaro 2/2
        Q = Q|set([k])
        #Marcamos a k como una tarea asignada
        j = int(sink)
        #Actualizamos las asignaciones 
        while i != k:
            i = pred[j]
            row[j] = i
            h = phi[i]
            phi[i] = j
            j = h

    for i in range(n):
        for j in range(n):
            #Obtenemos a la matriz de costos reducidos 
            #asociada a la solución dual
            C[i,j] = C[i,j]-u[i]-v[j]
            if row[j] == i:
            #Obtenemos la solución óptima 
                x[i,j] = 1
    #Obtenemos el valor óptimo
    z = sum(u)+sum(v)
    return(x,z)
    
def Value0(C,i,j):
    """
    Función auxiliar de BnB utilizada para ramificar sobre 
    x[i,j] = 0
    """
    C[i,j] = np.inf

def Value1(C,i,j):
    """
    Función auxiliar de BnB utilizada para ramificar sobre 
    x[i,j] = 1
    """
    n = len(C[0,:])
    for k in range(n):
        C[i,k] = np.inf
        C[k,j] = np.inf
    C[j,i] = np.inf
    C[i,j] = 0

def theta(C,i,j):
    """
    Función auxiliar de BnB, la cual regresa la penalización
    de cada índice [i,j] con respecto a una matriz reducida C.
    """
    a,b = np.inf,np.inf
    a = list(C[i,:])
    a.sort()
    b = list(C[:,j])
    b.sort()
    return a[1]+b[1]

def Branch(C):
    """
    Función auxiliar de BnB utilizada para obtener los 
    índices de ramificación [i,j] a partir de una matriz de 
    costos reducidos
    """
    Index = []
    n = len(C[0,:])
    for i in range(n):
        for j in range(n):
            if C[i,j] == 0:
            #Almacenamos a los índices [i,j] con costo
            #reducido 0
                Index.append([i,j])
    m = len(Index)
    values = m*[0]
    k = 0
    Hungarian(C)
    while k < m:
        values[k] = theta(C,Index[k][0],Index[k][1])
        if values[k] == np.inf:
        #Los descartamos si su penalización es infinito 
            values.pop(k)
            Index.pop(k)
            m += -1
        else:
            k += 1
    p = values.index(max(values))
    #Obtenemos el par [i,j] con mayor penalización finita
    return Index[p]
    
def constraint(X):
    """
    Función auxiliar de BnB utilizada para saber si X es 
    solución factible del PAV verificando que esta corresponda
    a un tour.
    """
    aux = 0        #El número de ciudades visitadas 
    index_i = 0    #El índice que nos indicará la fila
    index_j = 0    #El índice que nos indicará la columna
    n = len(X[0,:])#Número de ciudades
    Bool = True

    while index_i != 0 or Bool:
    #En este ciclo se verifica cuantas ciudades se visita con 
    #la matriz X a partir de la ciudad 0
        Bool = False
        while X[index_i, index_j] == 0:
            index_j += 1
        index_i = index_j
        index_j = 0
        aux += 1
    
    if aux == n:
    #Si resulta ser un recorrido de n ciudades, es factible
        return(True)
    else:
    #En otro caso, no lo es
        return(False)
        
def cap(lst1, lst2):
    """
    Convierte a las listas en conjuntos para encontrar la 
    lista intersectada (en el sentido de conjuntos). 
    """
    return list(set(lst1) & set(lst2))

def cup (lst1, lst2):
    """
    Convierte a las listas en conjuntos para encontrar la 
    lista unida (en el sentido de conjuntos).
    """
    return list(set(lst1)|set(lst2))

def Neighboor(node,A):
    """
    Dada una lista de arcos y un nodo, esta función regresa 
    un listado de los vecinos del nodo.
    """
    B = []
    m = len(A)
    for i in range(m):
        if A[i][0] == node:
            B.append(A[i][1])
        if A[i][1] == node:
            B.append(A[i][0])
    return(B)

def Cycles(n,A):
    """
    Suponiendo que una gráfica G1 = (X,A1) no contiene ciclos,
    este algoritmo nos ayuda a verificar que G2 = (X,A2) no 
    tenga, donde A2 = A1|{a} es la unión de dos conjuntos 
    ajenos. Esto lo hace al verificar que haya n-#A2 
    componentes conexas.
    """
    N = list(range(n)) #Listado de nodos de G
    E = [] #Listado auxiliar
    c = 0  #Contador de componentes conexas
    while len(N) != 0:
        #Mientras haya nodos sin examinar...
        i = N[0]    #Elegimos un nodo
        E.append(i) #Lo introducimos en la lista auxiliar
        while len(E) != 0:
            #Mientras el listado auxiliar tenga elementos
            i = E[0] #elegimos un elemento de E
            E.remove(i) 
            N.remove(i)
            #Lo retiramos de N y E, pues está siendo examinado
            E = cup(cap(N,Neighboor(i,A)),E)
            #Redefinimos a E con los vecinos de i
        c += 1 #Agregamos una componente conexa al contador
    if c == n - len(A):
    #Si hay n-#A componentes conexas, entonces no hay ciclos
        return False
    else:
    #En otro caso, A contiene ciclos
        return True
		
	

def NonFeasible(n, V1):
    """
    Función auxiliar de BnB, verifica que las variables 
    ramificadas con valor 1 no formen un circuito no 
    hamiltoniano, pues estos desembocan en nodos no factibles.
    """
    if len(V1) < n and Cycles(n,V1):
        return True
    else:
        return False


def BnB(C):
    """
    CUIDADO: ESTE ALGORITMO ES EXPONENCIAL, TOME 
    PRECAUCIONES DURANTE SU EJECUCIÓN CON INSTANCIAS GRANDES.
    Este algoritmo de ramificación y acotamiento resuelve
    instancias del problema del agente viajero codificadas en 
    una matriz de costos. Recordemos que el PAV consiste en 
    recorer un conjunto de n ciudades al menor costo pasando 
    sólo una vez por cada una y regresando a la ciudad 
    inicial (buscar el circuito hamiltoniano de menor costo).
    El algoritmo realiza una búsqueda exhaustiva implícita por
    utilizando como subrutina al algoritmo húngaro, el cual 
    resuelve instancias del PAV sin la restricción de 
    subtours, además ramifica a partir de un conjunto de 
    penalizaciones definido por Little, Murty, et al. en 1963.
    """
    i = 1     #índice de nodos totales
    index = 1 #índice del nodo con la mejor cota
    Ca = []   #Nodos candidatos
    NP = []   #Nodos no prometedores
    NF = []   #Nodos no factibles
    n = len(C[0,:])
    V0 = []   #índices de variables x's con valor 0
    V1 = []   #índices de variables x's con valor 1
    z0 = np.inf #El mejor valor encontrado hasta el momento
    X,w0 = Hungarian(C)
    #La solución relajada del PAV como problema de asignación
    #y su valor, el cual será una cota inferior del valor opt
    F = constraint(X)
    #F será verdadero iff X es un tour
    Nodes = {i:[w0, copy.deepcopy(X), copy.deepcopy(C), F, [np.inf,np.inf], V0, V1, [], "NT"]}
    #En nodes se almacenaran los nodos terminales y 
    Scout = np.array([[w0],[i]]) 
    #El arreglo Scout guardara los nodos (por índice) y sus
    #cotas para buscar por el nodo "más prometedor"
    while np.shape(Scout) != (2,0):
        index = Scout[1,0]
        #Consideramos al nodo con mejor cota
        print("La mejor cota actual es w = "+str(Scout[0,0]))
        if NonFeasible(n, Nodes[index][6]):
        #En caso de que las variables con valor 1 formen un
        #circuito no hamiltoniano, se etiqueta al nodo como
        #'No Factible'
            Nodes[index][8] = "NF"
            NF.append(index)
            Scout = np.delete(Scout,0,1)
        elif Nodes[index][0] <= z0 and Nodes[index][3]:
        #En caso de que el nodo tenga asociado un tour con 
        #al menos tan buen valor como z0...
            Nodes[index][8] = "Ca"
            #Lo etiquetamos como candidato
            z0 = Nodes[index][0]
            #Actualizamos el costo
            X = Nodes[index][1]
            #Guardamos el tour
            Ca.append(index)
            print("Hay "+str(len(Ca))+" nodos candidatos")
            print("El mejor valor encontrado es z = "+str(z0))
            Scout = np.delete(Scout,0,1)
        elif Nodes[index][0] > z0:
        #En caso de que el nodo tenga cota mayor al costo,
        #lo etiquetamos como "No Prometedor"
            Nodes[index][8] = "NP"
            NP.append(index)
            Scout = np.delete(Scout,0,1)
        else:
        #En caso de que el nodo no cumpla las condiciones
        #anteriores, se debe ramificar sobre él
            B = Branch(copy.deepcopy(Nodes[index][2]))
            #Se ramifica al índice [i,j] obtenido a partir de
            #las penalizaciones
            Nodes[index][4] = B
            ##################################################
            #Ramificación valuando a x[i,j] = 0
            C0 = copy.deepcopy(Nodes[index][2])
            Value0(C0,B[0],B[1])
            V0 = copy.deepcopy(Nodes[index][5])
            V0.append(B)
            X0,w0 = Hungarian(C0)
            print(w0+Nodes[index][0])
            F = constraint(X0)
            Nodes.update({(i+1):[w0+Nodes[index][0], copy.deepcopy(X0), copy.deepcopy(C0), F, B, V0, copy.deepcopy(Nodes[index][6]), index, "NT"]})
            ##################################################
            #Ramificación valuando a x[i,j] = 1
            C1 = copy.deepcopy(Nodes[index][2])
            Value1(C1,B[0],B[1])
            V1 = copy.deepcopy(Nodes[index][6])
            V1.append(B)
            X1,w1 = Hungarian(C1)
            F = constraint(X1)
            Nodes.update({(i+2):[w1+Nodes[index][0], copy.deepcopy(X1), copy.deepcopy(C1), F, B, copy.deepcopy(Nodes[index][5]), V1, index, "NT"]})
            ##################################################
            #Eliminamos al nodo actual de la lista Scout
            Scout = np.delete(Scout,0,1)
            NewNodes = np.array([[w0+Nodes[index][0],w1+Nodes[index][0]],[i+1,i+2]])
            #Agregamos a los nuevos nodos
            Scout = np.hstack((Scout,NewNodes.copy()))
            Scout = Scout[:, Scout[0].argsort()]
            [Nodes.pop(x) for x in [index]]
            i += 2
            #Actualizamos el contador de nodos
    return(X, z0)