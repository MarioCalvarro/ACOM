'''
GRUPO: Beñat Pérez de Arenaza Eizaguirre, Mario Calvarro Marines y Diego Ostos Arellano.
En la primera parte del problema, vamos a calcular la
multiplicación de matrices por medio del método de escuela. En este método
únicamente usamos la función multiplicación_mat_escuela(A,B,n), función a la
cual le pasamos las dos matrices A y B, así como la dimensión de la matriz. Esta
función es sencilla, vamos rellenando una matriz SOL por filas y columnas y en
cada posición (i,j) vamos añadiendo el valor de multiplicar los elementos de las
filas por los elementos de las columnas. De esta manera, obtenemos la matriz
solución de tamaño nXn. Para la segunda parte del problema, vamos a calcular la
multiplicación mediante el método de strassen para el cual nos basamos en la
descripción de pasos de link que es puesto en el enunciado. Por lo que siguiendo
estos pasos, realizamos lo siguiente. En primer lugar, si la dimensión de la
matriz es menor o igual a 5, se realiza la multiplicación de escuela, explicada
en el párrafo anterior. En caso contrario, continuamos con este método
mecionado. Si la matriz impar, tenemos que convertirla en par, lo que implica
que vamos a tener que añadir una columna y fila de 0 para hacerla así con un n
par. Seguidamente, procedemos a  “romper” las dos matrices en cuatro submatrices
(A11,A12,A21,A22 y B11,B12,B21,B22). Procedemos a realizar las sumas y restas
pertinentes como son descritas en el link asociado del ejercicio. Obteniendo así
las matrices: M1,M2,M3,M4,M5,M6 y M7. Estas 7 matrices son posteriormente usadas
para formar las cuatro submatrices (C11,C12,C21,C22) que tras una unión de
submatrices darán lugar a la matriz C. Pero está matriz todavía no es solución,
pues en el caso en el que al comienzo del procedimiento hayamos tenido que
añadir columna y fila de ceros, precisamos de eliminar esa fila y columna para
así obtener la matriz solución C que un matriz de dimensiones nXn como las
matrices A y B originales.
'''

import numpy as np
import random as rnd

def mul_mat_escuela(A, B, n):
    return [[sum([A[i][k] * B[k][j] for k in range(n)]) for j in range(n)] for i in range(n)]

def quitar_padding(matrix):
    n = len(matrix)
    return [row[:n-1] for row in matrix[:n-1]]

def suma_matrices(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def resta_matrices(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def mul_mat_strassen(A, B, n):
    if n <= 5:
        return mul_mat_escuela(A, B, n)

    # Si n es impar, añadimos una fila y columna de ceros
    padding = False
    if (n % 2 != 0):
        A = [r + [0] for r in A] + [[0] * (n + 1)]
        B = [r + [0] for r in B] + [[0] * (n + 1)]
        n += 1
        padding = True

    mid = n // 2

    # Separamos las matrices A y B en cuatro submatrices
    A11 = [row[:mid] for row in A[:mid]]
    A12 = [row[mid:] for row in A[:mid]]
    A21 = [row[:mid] for row in A[mid:]]
    A22 = [row[mid:] for row in A[mid:]]

    B11 = [row[:mid] for row in B[:mid]]
    B12 = [row[mid:] for row in B[:mid]]
    B21 = [row[:mid] for row in B[mid:]]
    B22 = [row[mid:] for row in B[mid:]]

    # El algoritmo de Strassen define estos valortes
    M1 = mul_mat_strassen(suma_matrices(A11, A22), suma_matrices(B11, B22), mid)
    M2 = mul_mat_strassen(suma_matrices(A21, A22), B11, mid)
    M3 = mul_mat_strassen(A11, resta_matrices(B12, B22), mid)
    M4 = mul_mat_strassen(A22, resta_matrices(B21, B11), mid)
    M5 = mul_mat_strassen(suma_matrices(A11, A12), B22, mid)
    M6 = mul_mat_strassen(resta_matrices(A21, A11), suma_matrices(B11, B12), mid)
    M7 = mul_mat_strassen(resta_matrices(A12, A22), suma_matrices(B21, B22), mid)

    # Combinación de los resultados
    C11 = suma_matrices(resta_matrices(suma_matrices(M1, M4), M5), M7)
    C12 = suma_matrices(M3, M5)
    C21 = suma_matrices(M2, M4)
    C22 = suma_matrices(resta_matrices(suma_matrices(M1, M3), M2), M6)

    # Unimos las submatrices
    C = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(mid):
        for j in range(mid):
            C[i][j] = C11[i][j]
            C[i][j + mid] = C12[i][j]
            C[i + mid][j] = C21[i][j]
            C[i + mid][j + mid] = C22[i][j]

    if padding:
        return quitar_padding(C)
    else:
        return C


def compare_with_np(A, B):
    A_np = np.array(A)
    B_np = np.array(B)

    C_np = np.matmul(A_np, B_np)

    C_pre = mul_mat_strassen(A, B, len(A))

    C = np.array(C_pre)
    return (C==C_np).all()

def generar_matrices(m):
    A = [[float(rnd.randint(1, 100)) for _ in range(m)] for _ in range(m)]

    B = [[float(rnd.randint(1, 100)) for _ in range(m)] for _ in range(m)]
    return A, B

def probar(m):
    A, B = generar_matrices(m)
    return compare_with_np(A, B)
