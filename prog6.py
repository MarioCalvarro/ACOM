# Mario Calvarro, Beñat Pérez
def rotar(f, t):
    if t == 0:
        return f

    l = len(f)
    res = [0]*l

    for j in range(l):
        if j + t >= l or j + t < 0:
            # [u^n2] = [-1]
            res[( j + t ) % l] = -f[j]
        else:
            res[( j + t ) % l] = f[j]

    return res

def sumar_vectores(f, g, p):
    h = []
    for i in range(len(f)):
        h += [(f[i] + g[i]) % p]
    return h

def restar_vectores(f, g, p):
    h = []
    for i in range(len(f)):
        h += [(f[i] - g[i]) % p]
    return h

# Algoritmo de Cooley-Tuckey para calcular DFT_n(p)
# n = len(p) debe ser una potencia de 2 y xi debe ser
# una raíz n-ésima primitiva de la unidad
def fft_adaptado(f, xi, p):
    n = len(f)
    if n == 1:
        return f
    p_even = [0j] * (n//2)
    p_odd  = [0j] * (n//2)
    for i in range(n//2):
        p_even[i] = f[2*i]
        p_odd[i]  = f[2*i+1]

    a_even = fft_adaptado(p_even, xi*2, p)
    a_odd  = fft_adaptado(p_odd,  xi*2, p)

    a = [0] * n
    for i in range(n//2):
        rotar_impar = rotar(a_odd[i], i*xi)

        a[i] = sumar_vectores(a_even[i], rotar_impar, p)
        a[i + n // 2] = restar_vectores(a_even[i], rotar_impar, p)
    return a

# Esta función calcula la transformada inversa
# utilizando que D(xi)^(-1) = 1/n * D(1/xi)
def ifft_adaptado(a, xi, p):
    n = len(a)
    transformado = fft_adaptado(a, 1.0/xi, p)
    for i in range(n):
        transformado[i] /= n
    return transformado


def split(f, k):
    """Divide un polinomio en representación con u y v."""
    n = len(f)
    k2 = (k + 1) // 2 # 2
    k1 = k - k2       # 1
    n1 = 2**k1
    #TODO: Asegurar que los ceros se añaden bien
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

def negaconvolucion(f, g, p):
    return [0]

def mult_ss_mod(f, g, k, p):
    if f == [] or g == []:
        return []

    # Casos base:
    if k == 0:
        return [(f[0]*g[0]) % p]

    if k == 1:
        return [(f[0]*g[0] - f[1]*g[1])%p, (f[0]*g[1] + f[1]*g[0]) % p]

    if k == 2:
        a_0 = (f[0]*g[0] - f[1]*g[3] - f[2]*g[2] - f[3]*g[1]) % p
        a_1 = (f[0]*g[1] + f[1]*g[0] - f[2]*g[3] - f[3]*g[2]) % p
        a_2 = (f[0]*g[2] + f[1]*g[1] + f[2]*g[0] - f[3]*g[3]) % p
        a_3 = (f[0]*g[3] + f[1]*g[2] + f[2]*g[1] + f[3]*g[0]) % p
        return [a_0, a_1, a_2, a_3]

    f_class = split(f, k)
    g_class = split(g, k)

    h_class = negaconvolucion(f_class, g_class, p)
    return join(h_class)

def eliminar_ceros(v):
    n = len(v)
    while n >= 1 and v[n-1] == 0:
        n -= 1
    del v[n:]
    return v


def mult_pol_mod(f, g,p):
    if f == [] or g == []:
        return []

    m = len(f) + len(g)
    k = 0
    while m != 0:
        m //= 2
        k += 1

    # Añadimos los ceros necesarios
    f = f + [0]*(2**k - len(f))
    g = g + [0]*(2**k - len(g))

    h = eliminar_ceros(mult_ss_mod(f, g, k, p))
    return h
