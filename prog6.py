# Mario Calvarro, Beñat Pérez
def rotar(f, t):
    if t == 0:
        return f

    n = len(f)
    res = [0] * n

    for i in range(n):
        if i + t >= n or i + t < 0:
            # [u^n2] = [-1]
            res[( i + t ) % n] = -f[i]
        else:
            res[( i + t ) % n] = f[i]

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

    p_even = [0] * (n//2)
    p_odd  = [0] * (n//2)
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
    n1 = len(a[0])
    n2 = len(a)
    transformado = fft_adaptado(a, -xi, p)
    inversa = pow(n2, -1, p)

    for i in range(n2):
        for j in range(n1):
            transformado[i][j] = (transformado[i][j] * inversa) % p

    return transformado


def split(f, k):
    """Divide un polinomio en representación con u y v."""
    k2 = (k + 1) // 2
    k1 = k - k2

    n1 = 2**k1
    n2 = 2**k2
    res = []
    for i in range(0, n2):
        res += [f[i*n1 : i*n1 + n1] + [0]*n1]

    return res

def join(f_split, p):
    """Reconstruye un polinomio desde su representación con u y v."""
    n2 = len(f_split)
    n1 = len(f_split[0])//2

    pol = restar_vectores(f_split[0][0:n1], f_split[n2-1][n1:(2*n1)], p)
    for i in range(0,n2-1):
        pol   = pol + sumar_vectores(f_split[i][n1:(2*n1)], f_split[i+1][0:n1], p)

    return pol
    # m = len(f_split[0])
    # k = len(f_split) * len(f_split[0])
    # f = [0] * k
    # for i in range(len(f_split)):
    #     for j in range(len(f_split[i])):
    #         f[i * m + j] = f_split[i][j]
    # return f

def negaproducto(f, n1, n2):
    """Realiza el producto w·f de la negaconvolucion. En este caso, son
    simplemente rotaciones de los elementos de f"""
    t = (2*n1) // n2
    res = [f[0]]

    for i in range (1, n2):
        res += [rotar(f[i], i*t)]

    return res

def inv_negaproducto(f, n1, n2):
    """Realiza el producto w·f de la negaconvolucion. En este caso, son
    simplemente rotaciones inversas de los elementos de f"""
    t = (2*n1) // n2
    res = [f[0]]

    for i in range (1, n2):
        res += [rotar(f[i], -(i*t))]

    return res

def negaconvolucion(f, g, k1, p):
    n1 = len(f[0])//2       #TODO: Bien?
    n2 = len(f)

    term1 = fft_adaptado(negaproducto(f, n1, n2), (4*n1)//n2, p)
    term2 = fft_adaptado(negaproducto(g, n1, n2), (4*n1)//n2, p)

    prod = []
    for i in range(n2):
        prod += [mul_ss_pol_mod(term1[i], term2[i], k1 + 1, p)]

    transformado_inv = ifft_adaptado(prod, (4*n1)//n2, p)

    return inv_negaproducto(transformado_inv, n1, n2)

def mul_ss_pol_mod(f, g, k, p):
    if f == [] or g == []:
        return []

    #Casos base:
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


    k2 = (k + 1) // 2
    k1 = k - k2

    f_class = split(f, k)
    g_class = split(g, k)

    h_class = negaconvolucion(f_class, g_class, k1, p)

    return join(h_class, p)

def eliminar_ceros(v):
    n = len(v)
    while n >= 1 and v[n-1] == 0:
        n -= 1
    del v[n:]
    return v


def mul_pol_mod(f, g, p):
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

    h = eliminar_ceros(mul_ss_pol_mod(f, g, k, p))
    return h

