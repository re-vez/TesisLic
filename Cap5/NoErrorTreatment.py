#José Fernando Méndez Torres
#Códigos correspondientes al quinto capítulo de la tesis  
#"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.

#Este código realiza el tratamiento de datos descrito por 
#Gusfield en Integer Linear Programming in Computational and
#System Biology, subsección 9.3.2 para obtener la matriz de 
#pesos a partir de los traslapes de lecturas

Length = 0
#Esta variable nos ayudará a calcular la cobertura

Count = 1
#El contador nos ayudará a identificar las líneas donde hay 
#lecturas
Fragments = []
#Aquí se almacenarán las lecturas
with open("J02459.fastq","r") as f:
    for line in f:
        if Count%4 == 2:
        #Recordemos que si el número de la línea módulo 4 es 
        #igual a 2, hay una lectura
            line = line.replace('\n', '')
            #Eliminamos el salto de línea de la cadena
            Fragments.append(line)
            #La almacenamos
        Count += 1

##############################################################
#Eliminación de lecturas cortas
##############################################################
#Fijamos al índice i en 0
i = 0
#Mientras haya fragmentos sin revisar...
while i < len(Fragments):
    l = len(Fragments[i])
    if l < 100:
    #Si la longitud del fragmento es menor a 100
        Fragments.pop(i)
        #Lo desechamos
    else:
        #En otro caso, 
        Length += l
        #Sumamos la longutud a Length
        i += 1
        #Pasamos al siguiente índice

##############################################################
#Eliminación de lecturas contenidas en otras
##############################################################
#Fijamos a los índices auxiliares i y j para comparar todas 
#las lecturas de Fragments cuidando que j < i
i = 1, j = 0
while i < len(Fragments):
    while j < i and i < len(Fragments):
    #Mientras no se hayan revisado los pares (i,j) con i < j
        FragmentI = Fragments[i]
        FragmentJ = Fragments[j]
        CondJinI = FragmentI.find(FragmentJ)
        CondIinJ = FragmentJ.find(FragmentI)
        #Verificamos que la lectura i no esté dentro de la 
        #lectura j y viceversa
        if CondJinI >= 0:
        #Si la j-ésima lectura resulta estar contenida en la 
        #i-ésima
            Length -= len(FragmentJ)
            #Dejamos de considerarla en Length
            Fragments.pop(j)
            #La eliminamos
            i -= 1
            #Retrocedemos un índice i
        elif CondIinJ >= 0:
        #Si la i-ésima lectura resulta estar contenida en la 
        #j-ésima
            Length -= len(FragmentI)
            #Dejamos de considerarla en Length
            Fragments.pop(i)
            #La eliminamos
        else:
            #Avanzamos con el índice j
            j+=1
    #Restauramos a j = 0 y avanzamos sobre i una unidad
    j = 0
    i += 1
#Imprimimos la cantidad de bases que secuenciamos 
print(Length)

#Guardamos las lecturas tratadas en un archivo .txt
with open("Treated.txt","w") as f:
    for i in range(len(Fragments)):
        f.write(Fragments[i]+"\n")
