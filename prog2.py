'''
Para la resolución de este ejercicio hemos usado programación dinámica.
En primer lugar hemos inicializado una matriz de nxmx3 a una cota superior de los posibles valores que
podría tomar una cierta casilla. Para ello hemos sumado todos los valores positivos de la tabla.
Para cada casilla tenemos 3 un array de 3 elementos, cada uno de los cuales refleja
el valor mínimo al que se puede llegar a dicha casilla habiendo hecho uso de 0, 1 o 2 movimientos
consecutivos de caballo como últimos movientos respectivamente. Es decir, la mejor forma de llegar
a la casilla [2,2] sin que el moviento previo sea un movimiento de caballo es 5, siendo el último 
movimiento de caballo pero el inmediatamente anterior no sea de caballo es 6, y si los dos movimientos 
anteriores sean de caballo y el anterior a estos dos no es 7, entonces el array matriz[2][2]
sería [5,6,7].

Sabiendo esto, para la resolución del problema, hemos hecho un bucle que recorra todas las casillas 
menos la [0][0] que hemos inicializado al valor de a[0][0], que será nuestro caso base. Esto es así porque
todos los posibles itinerarios desde [0][0] hasta [n-1][m-1] pasan por ambas dos casillas. 
En el bucle, siempre respetando que las casillas que se van a consultar existan, vamos a calcular
el valor mínimo que se necesita para llegar a la casilla [i,j] apoyándonos en las 3 casillas de donde
puede estar viniendo la ficha. Así, sumamos el valor de la casilla actual al mínimo de las tres posibles
casillas anteriores desde las que ha podido saltar a la actual. Todo ello lo hacemos para los 3 casos
explicados anteriormente (haber llegado a la casilla actual sin un salto de caballo, llegar con un salto,
y llegar con dos saltos consecutivos).
Finalmente, devolvemos el valor mínimo de llegar a la casilla que está abajo a la derecha, para lo cual
calculamos el mínimio entre las 3 formas distintas (desde el punto de vista del número de movimientos
de caballo consecutivos que se han hecho inmediatamente antes de llegar a la casilla actual)
de haber llegado a dicha casilla.
'''

def min_suma_casillas(n,m,a):
    maximo = 0
    for i in range(n): # calculamos una cota superior
        for j in range(m):
            if a[i][j] > 0:
                maximo = maximo + a[i][j]

    matriz = [[[maximo] * 3 for _ in range(m)] for _ in range(n)]    
    matriz[0][0][0] = a[0][0]

    for i in range(n):
        for j in range(m):

            if i == 0 and j == 0:
                continue

            arriba, izquierda = maximo, maximo
            
            if i > 0: #viene de arriba
                arriba = min(matriz[i-1][j])

            if j > 0: # viene de la izquierda
                izquierda = min(matriz[i][j-1])

            if j > 1 and i > 0: #hace un salto de caballo
                matriz[i][j][1] = matriz[i-1][j-2][0] + a[i][j]
                matriz[i][j][2] = matriz[i-1][j-2][1] + a[i][j]

            matriz[i][j][0] = min(arriba, izquierda) + a[i][j]
           

    return min (matriz[n-1][m-1])
