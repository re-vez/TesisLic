# Introducción

El ADN es la molécula que contiene las instrucciones para generar proteínas, las cuales se encargan del crecimiento, desarrollo, funcionalidad y reproducción de los seres vivos y es representada mediante una **secuencia** o sucesión de cuatro letras, donde cada letra corresponde a un compuesto molde llamado base nitrogenada [1]. Se denomina **secuenciación** a la determinación del orden de las cuatro diferentes bases de alguna de estas moléculas, es decir, la determinación de la secuencia [2]. Para ello, existen múltiples técnicas que han sido desarrolladas a lo largo de las últimas décadas. La cuestión es que dichos métodos sólo pueden obtener secuencias cortas de ADN. Conocer el orden de las bases en estos compuestos nos permite saber cómo funcionan los organismos a nivel molecular, por lo que es de gran interés para las ciencias biológicas y de la salud.

Si poseemos la secuencia de algún individuo de cierta especie, obtener la de cualquier otro miembro resulta más fácil, ya que se tiene una referencia para las demás. Como menciona Ben Langmead [3], los métodos de secuenciación nos brindan las piezas de un rompecabezas y, por ende, la molécula de referencia funge como la imagen en la caja que nos indica cómo armarlo. El problema surge en cómo obtener dicha secuencia, es decir ¿cómo obtenemos la imagen de la caja? Una técnica para obtenerla es a través del proceso de **secuenciación por escopeta**, el cual consiste de un procedimiento experimental y un análisis informático [4, pág. 159].

En los siglo XIX y XX, era común que los vendedores se transportaran a varias ciudades. Dado un listado de ciudades a recorrer y un costo de viaje entre cada ciudad, buscaban trasladarse por todas ellas al menor costo posible, pasando sólo una vez a cada una y regresando a su ciudad de origen [5]. El problema anterior es conocido como el **problema del agente viajero (PAV)** y ha sido altamente estudiado desde la década de 1930 hasta la fecha.

Aunque de primera vista el problema del agente viajero y la secuenciación por escopeta parecieran no presentar relación alguna, se puede utilizar el planteamiento del primero en el análisis informático del segundo para así obtener la secuencia de referencia, o la antes mencionada imagen de la caja. El objetivo principal de este trabajo es estudiar el ordenamiento de secuencias de ADN modelándolo matemáticamente como un problema del agente viajero. Este texto está dividido en 5 capítulos donde se tratan temas de biología, matemáticas y computación.

En el primer capítulo veremos algunos conceptos básicos de biología molecular para poder describir a las proteínas y al ADN. También veremos el dogma central de la biología molecular, el cual es una simplificación del flujo de información de cómo se sintetizan las proteínas a partir del ADN dentro de las células. En este mismo capítulo se presentan dos métodos utilizados para conocer la secuencia de pequeñas moléculas de ADN.

En el segundo capítulo plantearemos el problema del agente viajero de dos formas, como una red y como un modelo de programación entera mixto. Veremos el problema de asignación, el cual es una versión relajada del problema entero mixto y estudiaremos un algoritmo que lo resuelve eficientemente, llamado algoritmo húngaro. Además, proporcionaremos un algoritmo de ramificación y acotamiento que resuelve el PAV.

En el tercer capítulo examinaremos los errores comunes de secuenciación, el manejo y tratamiento de los archivos obtenidos de las máquinas secuenciadoras y el planteamiento del PAV a partir de dichos datos, con el cual obtendremos la molécula original de ADN. El modelo es desarrollado de dos formas, primero suponiendo que no tenemos errores en la parte experimental y posteriormente considerando dichas fallas.

En el cuarto capítulo daremos la definición de algoritmo a partir del modelo computacional de Turing, de ahí justificaremos que el problema del agente viajero es difícil de resolver desde la teoría de la complejidad, es decir, demostraremos que es **NP-Duro**. Asimismo, introduciremos brevemente a los algoritmos heurísticos y desarrollaremos un algoritmo genético para el PAV con el que se pretende encontrar buenas soluciones.

En el quinto capítulo trabajaremos con el ADN de referencia del fago lambda. A partir de él, generaremos datos sin errores por medio de una computadora, aplicaremos el modelo propuesto en el capítulo 3 y lo secuenciaremos con los algoritmos presentados en los capítulos 2 y 4 para comparar los resultados con la secuencia original.

Referencias

1. H. Curtis, H. S. Barnes, A. Schnek y A. Massarini, Biología, 7.a ed., E. M. Panamericana, ed. 2008, págs. 172-204.
2. J. M. Heather y B. Chain, «The sequence of sequencers: The history of sequencing DNA», Genomics, vol. 107, n.o 1, págs. 1-8, 2016. doi:10.1016/j.ygeno.2015.11.003.
3. B. Langmead. (2015). Sequencers give pieces to genomic puzzles, Youtube, dirección: https://youtu.be/X307VyAeHfI.
4. D. Gusfield, Integer Linear Programming in Computational and Systems Biology. Cambridge University Press, 2019.
5. D. L. Applegate, R. E. Bixby, V. Chvatál y W. J. Cook, The Traveling Salesman Problem: A Computational Study. Princeton University Press, 2006, págs. 1-4.


