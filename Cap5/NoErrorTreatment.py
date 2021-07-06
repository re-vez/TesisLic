#José Fernando Méndez Torres
#Códigos correspondientes al quinto capítulo de la tesis  
#"Secuenciación del ADN por medio del problema del agente 
#viajero". Los métodos son descritos a detalle en dicha tesis.

#Este código realiza el tratamiento de datos descrito por 
#Gusfield en Integer Linear Programming in Computational and
#System Biology, subsección 9.3.2 para obtener la matriz de 
#pesos a partir de los traslapes de lecturas
Length = 0

Count = 1
Fragments = []
with open("J02459.fastq","r") as f:
    for line in f:
        if Count%4 == 2:
            line = line.replace('\n', '')
            Fragments.append(line)
        Count += 1


i = 0

while i < len(Fragments):
    l = len(Fragments[i])
    if l < 100:
        Fragments.pop(i)
    else:
        Length += l
        i += 1

i = 1
j = 0

while i < len(Fragments):
    while j < i and i < len(Fragments):
        FragmentI = Fragments[i]
        FragmentJ = Fragments[j]
        CondJinI = FragmentI.find(FragmentJ)
        CondIinJ = FragmentJ.find(FragmentI)
        if CondJinI >= 0:
            Fragments.pop(j)
            i -= 1
        elif CondIinJ >= 0:
            Fragments.pop(i)
        else:
            j+=1
    j = 0
    i += 1
print(Length)

with open("Treated.txt","w") as f:
    for i in range(len(Fragments)):
        f.write(Fragments[i]+"\n")
