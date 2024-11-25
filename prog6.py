# Mario Calvarro, Beñat Pérez
# Algoritmo de Cooley-Tuckey para calcular DFT_n(p)
# n = len(p) debe ser una potencia de 2 y xi debe ser
# una raíz n-ésima primitiva de la unidad
def fft(p, xi):
    n = len(p)
    if n == 1:
        return p
    p_even = [0j] * (n//2)
    p_odd  = [0j] * (n//2)
    for i in range(n//2):
        p_even[i] = p[2*i]
        p_odd[i]  = p[2*i+1]
    a_even = fft(p_even, xi**2)
    a_odd  = fft(p_odd,  xi**2)
    a = [0j] * n
    for i in range(n//2):
        a[i]      = a_even[i] + xi**i * a_odd[i]
        a[i+n//2] = a_even[i] - xi**i * a_odd[i]
    return a

# Esta función calcula la transformada inversa
# utilizando que D(xi)^(-1) = 1/n * D(1/xi)
def ifft(a, xi):
    n = len(a)
    p = fft(a, 1.0/xi)
    for i in range(n):
        p[i] /= n
    return p


def split(f, k):
    """Divide un polinomio en representación con u y v."""
    n = len(f)
    k2 = (k + 1) // 2 # 2
    k1 = k - k2       # 1
    n1 = 2**k1
    return [[f[i] for i in range(j, j + n1)] for j in range(0, n, n1)]

def join(f_split):
    """Reconstruye un polinomio desde su representación con u y v."""
    m = len(f_split[0])
    k = len(f_split) * len(f_split[0])
    f = [0] * k
    for i in range(len(f_split)):
        for j in range(len(f_split[i])):
            f[i * m + j] = f_split[i][j]
    return f

