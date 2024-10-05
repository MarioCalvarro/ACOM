'''
GRUPO: Beñat Pérez de Arenaza Eizaguirre, Mario Calvarro Marines y Diego Ostos Arellano
En la primera parte del problema, vamos a calcular la multiplicación de matrices por medio del método de escuela. En este método únicamente usamos 
la función multiplicación_mat_escuela(A,B,n), función a la cual le pasamos las dos matrices A y B, así como la dimensión de la matriz. Esta función
es sencilla, vamos rellenando una matriz SOL por filas y columnas y en cada posición (i,j) vamos añadiendo el valor de multiplicar los elementos de
las filas por los elementos de las columnas. De esta manera, obtenemos la matriz solución de tamaño nXn.
Para la segunda parte del problema, vamos a calcular la multiplicación mediante el método de strassen para el cual nos basamos en la descripción de
pasos de link que es puesto en el enunciado. Por lo que siguiendo estos pasos, realizamos lo siguiente. En primer lugar, si la dimensión de la 
matriz es menor o igual a 5, se realiza la multiplicación de escuela, explicada en el párrafo anterior. En caso contrario, continuamos con este
método mecionado. Si la matriz impar, tenemos que convertirla en par, lo que implica que vamos a tener que añadir una columna y fila de 0 para
hacerla así con un n par. Seguidamente, procedemos a  “romper” las dos matrices en cuatro submatrices (A11,A12,A21,A22 y B11,B12,B21,B22).
Procedemos a realizar las sumas y restas pertinentes como son descritas en el link asociado del ejercicio. Obteniendo así las
matrices: M1,M2,M3,M4,M5,M6 y M7. Estas 7 matrices son posteriormente usadas para formar las cuatro submatrices (C11,C12,C21,C22) que tras una
unión de submatrices darán lugar a la matriz C. Pero está matriz todavía no es solución, pues en el caso en el que al comienzo del procedimiento
hayamos tenido que añadir columna y fila de ceros, precisamos de eliminar esa fila y columna para así obtener la matriz solución C que un 
matriz de dimensiones nXn como las matrices A y B originales. 
'''
def mul_mat_escuela(A,B,n):
    SOL = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                SOL[i][j] += A[i][k]*B[k][j]
    return SOL

def sumar(A,B):
    n = len(A)
    return [[A[i][j] + B[i][j] for j in range(n)] for i in range(n)]

def restar(A,B):
    n = len(A)
    return [[A[i][j] - B[i][j] for j in range(n)] for i in range(n)]

def seperar_matriz(A,n):
    centro = n //2
    A11 = [row[:centro] for row in A[:centro]]
    A12 = [row[centro:] for row in A[:centro]]
    A21 = [row[:centro] for row in A[centro:]]
    A22 = [row[centro:] for row in A[centro:]]
    return A11,A12,A21,A22

def unir_matrices(C11,C12,C21,C22,n):
    SOL = [[0 for _ in range(n)] for _ in range(n)]
    centro = n//2
    for i in range(n):
        for j in range(n):
            if i < centro: 
                if j < centro: 
                   SOL[i][j] = C11[i][j]
                else:
                   SOL[i][j] = C12[i][j - centro]
            else:
                if j < centro: 
                   SOL[i][j] = C21[i-centro][j]
                else:
                   SOL[i][j] = C22[i-centro][j - centro]    
    return SOL

           


def mul_mat_strassen(A, B, n):
    if n <= 5:
        return mul_mat_escuela(A, B, n)
    
    # Si n es impar, añadimos una fila y columna de ceros
    impar = False
    if n % 2 != 0:
        A = [row + [0] for row in A] + [[0] * (n + 1)]
        B = [row + [0] for row in B] + [[0] * (n + 1)]
        n += 1
        impar = True

    # Separamos las matrices A y B en cuatro submatrices
    A11, A12, A21, A22 = seperar_matriz(A,n)
    B11, B12, B21, B22 = seperar_matriz(B,n)

    # El algoritmo de Strassen define estos valortes
    M1 = mul_mat_strassen(sumar(A11, A22), sumar(B11, B22), n // 2)
    M2 = mul_mat_strassen(sumar(A21, A22), B11, n // 2)
    M3 = mul_mat_strassen(A11, restar(B12, B22), n // 2)
    M4 = mul_mat_strassen(A22, restar(B21, B11), n // 2)
    M5  = mul_mat_strassen(sumar(A11, A12), B22, n // 2)
    M6 = mul_mat_strassen(restar(A21, A11), sumar(B11, B12), n // 2)
    M7 = mul_mat_strassen(restar(A12, A22), sumar(B21, B22), n // 2)

    # Combinación de los resultados
    C11 = sumar(restar(sumar(M1, M4), M5), M7)
    C12 = sumar(M3, M5 )
    C21 = sumar(M2, M4)
    C22 = sumar(restar(sumar(M1, M3), M2), M6)

    # Unimos las submatrices
    C = unir_matrices(C11, C12, C21, C22,n)

    # Remover la fila y columna extra si n era impar
    if impar:
        C = [row[:-1] for row in C[:-1]]
   

    return C