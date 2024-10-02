import numpy as np
import random as rnd

def mul_mat_escuela(A, B, n):
    return [[sum([A[i][k] * B[k][j] for k in range(n)]) for j in range(n)] for i in range(n)]


def add_padding(matrix):
    n = len(matrix)
    for i in range(n):
        matrix[i].append(0);
    matrix.append([0 for _ in range(n + 1)])


def remove_padding(matrix):
    n = len(matrix)
    return [row[:n-1] for row in matrix[:n-1]]

def suma_matrices(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def resta_matrices(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def mul_mat_strassen(A, B, n):
    """Perform matrix multiplication using Strassen's algorithm."""
    if n <= 5:
        return mul_mat_escuela(A, B, n)

    padding = False
    if (n % 2 != 0):
        # Make matrices even-sized if necessary
        # add_padding(A)
        # add_padding(B)
        A = [r + [0] for r in A] + [[0] * (n + 1)]
        B = [r + [0] for r in B] + [[0] * (n + 1)]
        n += 1
        padding = True

    mid = n // 2

    # Splitting the matrices into quadrants
    A11 = [row[:mid] for row in A[:mid]]
    A12 = [row[mid:] for row in A[:mid]]
    A21 = [row[:mid] for row in A[mid:]]
    A22 = [row[mid:] for row in A[mid:]]

    B11 = [row[:mid] for row in B[:mid]]
    B12 = [row[mid:] for row in B[:mid]]
    B21 = [row[:mid] for row in B[mid:]]
    B22 = [row[mid:] for row in B[mid:]]

    # Strassen's algorithm
    M1 = mul_mat_strassen(suma_matrices(A11, A22), suma_matrices(B11, B22), mid)
    M2 = mul_mat_strassen(suma_matrices(A21, A22), B11, mid)
    M3 = mul_mat_strassen(A11, resta_matrices(B12, B22), mid)
    M4 = mul_mat_strassen(A22, resta_matrices(B21, B11), mid)
    M5 = mul_mat_strassen(suma_matrices(A11, A12), B22, mid)
    M6 = mul_mat_strassen(resta_matrices(A21, A11), suma_matrices(B11, B12), mid)
    M7 = mul_mat_strassen(resta_matrices(A12, A22), suma_matrices(B21, B22), mid)

    # Combining the results into a single matrix
    C11 = suma_matrices(resta_matrices(suma_matrices(M1, M4), M5), M7)
    C12 = suma_matrices(M3, M5)
    C21 = suma_matrices(M2, M4)
    C22 = suma_matrices(resta_matrices(suma_matrices(M1, M3), M2), M6)

    # Recombining the four quadrants into one result matrix
    C = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(mid):
        for j in range(mid):
            C[i][j] = C11[i][j]
            C[i][j + mid] = C12[i][j]
            C[i + mid][j] = C21[i][j]
            C[i + mid][j + mid] = C22[i][j]

    if padding:
        return remove_padding(C)
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
    A = [[rnd.randint(0, 100) for _ in range(m)] for _ in range(m)]

    B = [[rnd.randint(0, 100) for _ in range(m)] for _ in range(m)]
    return A, B

def probar(m):
    A, B = generar_matrices(m)
    return compare_with_np(A, B)
