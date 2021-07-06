#José Fernando Méndez Torres
#Códigos correspondientes al quinto capítulo de la tesis  
#"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.

#En este archivo se simula la secuenciación de lecturas de un 
#genoma. El código genera un archivo .fastq a partir de un 
#archivo .fasta o .txt. 
import numpy as np

# Fijamos la semilla del proceso pseudo-aleatorio
np.random.seed(13013)

# Manipulación de parámetros 
PCRCycles = 5 #Para C ciclos PCR, habrá 2**C moléculas de ADN
LecturesLength = 200 #La longitud de las lecturas
Quality = 40 #Este parámetro no es conciderado en el proceso
             #La calidad de las lecturas
MeanDNaseCut = 200 #La longitud promedio de los cortes en las 
                   #Moléculas de ADN

Q = chr(Quality+64)
p = 1/(1+MeanDNaseCut)

#Datos del experimento in silico
SpeciesName = "Escherichia_phage_Lambda"
Code = "J02459"

#Obtención de la cadena de ADN
with open("Escherichia_phage_Lambda.txt","r") as f:
    _ = f.readline()
    DNA = '' #Inicializamos una cadena vacía
    for line in f:
    #Para cada linea en el archivo
        line = line.replace('\n', '')
        #Retiramos el salto de línea
        DNA += line
        #Guardamos la línea al final de DNA
#Obtenemos la longitud del genoma
GenomeLength = len(DNA)

#Generación de sitios de corte pseudo-aleatorios
GeometricNums = []
#Inicializamos esta lista para guardar los n(2**C) cortes 
#aleatorios de las 2**C moléculas
for i in range(2**PCRCycles):
#Para cada molécula
    GenomeCuts = [0] 
    #Aquí almancenaremos los cortes de la i-ésima molécula
    DNaseCutsAux = 0 
    #Esta variable auxiliar nos ayudará a mantener el 
    #seguimiento de los saltos de la función ( k = n(i) )
    while DNaseCutsAux < GenomeLength:
    #Mientras no hayamos recorrido todo el genoma con los #números aleatorios
        RandomCut = np.random.geometric(p)+DNaseCutsAux
        #Generamos un nuevo corte aleatorio
        if RandomCut > GenomeLength:
        #Si el corte traspasó la longitud del genoma
            DNaseCutsAux = GenomeLength
            #Lo hemos recorrido todo
        else:
            DNaseCutsAux = RandomCut
            #En otro caso, guardamos el sitio de corte en la 
            #variable auxiliar 
        GenomeCuts.append(DNaseCutsAux)
        #Guardamos el número aleatorio
    GeometricNums.append(GenomeCuts)
    #Guardamos los cortes de la molécula



#Generador de lecturas
with open("J02459.fastq","w") as f:
    Count = 1
    #Este contador asignará un índice a las lecturas generadas
    for i in range(2**PCRCycles):
    #Para todas las moléculas i...
        for j in range(len(GeometricNums[i])-1):
        #Para todos los cortes aleatorios j...
            Begin = GeometricNums[i][j]
            #Obtenemos el índice donde empieza la lectura
            n = GeometricNums[i][j+1]-GeometricNums[i][j]
            #Y la distancia entre el corte actual y su sucesor
            Length = min(LecturesLength, n)
            #Recordemos que simulamos LecturesLength ciclos 
            #de la SBS, por lo que la longitud de la lectura
            #está dada por el mínimo de esos dos números
            f.write("@"+Code+"."+str(Count)+" "+SpeciesName+" length="+str(Length)+"\n")
            #Escribimos la información de la lectura y el 
            #experimento
            f.write(DNA[Begin:Begin+Length]+"\n")
            #Plasmamos la lectura
            f.write("+"+Code+"."+str(Count)+" "+SpeciesName+" length="+str(Length)+"\n")
            #Reescribimos la información anterior
            f.write(Q*Length+"\n")
            #Plasmamos la calidad de la lectura
            Count += 1
