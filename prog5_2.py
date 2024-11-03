import random

# ---------------- FUNCIONES PREVIAS QUE USAREMOS ----------------

def jacobi(a, N):
    """Devuelve el símbolo de Jacobi usando recursión de cola"""
    m = 1
    while a != 0 and a != 1:
        if a < 0:
            if N % 4 == 1:
                m *= -1
            (a, N) = (-a, N)
        elif a >= N:
            (a, N) = (a % N, N)
        elif a % 2 == 0:
            if N % 8 == 1 or N % 8 == 7:
                (a, N) = (a // 2, N)
            else:
                m *= -1
                (a, N) = (a // 2, N)
        else:
            if N % 4 == 1 or a % 4 == 1:
                (a, N) = (N, a)
            else:
                m *= -1
                (a, N) = (N, a)
    if a == 0:
        m *= 0
    return m

def desc(p):
    """Para p primo, devuelve tupla (p, r, k) tal que p = r * 2**k + 1"""
    r = p - 1
    k = 0
    while r % 2 == 0:
        k += 1
        r //= 2
    return r, k

def orden(a, n):
    """Calcula el orden del elemento a en un grupo de n elementos"""
    for i in range(1, n):
        if pow(a, i, n) == 1:
            return i

def gen_raiz_prim(p):
    """Genera una raíz primitiva aleatoria para p primo"""
    a = random.randint(1, p)
    while orden(a, p) != p - 1:
        a = random.randint(1, p)
    return a

# ---------------- CÓDIGO PARA RESOLVER x^2 = a (p**N) ----------------
# --------- Resolvemos primero caso n == 1

def sqrt_mod1(a, p):
    # por cómo se plantea el problema puede que a no esté módulo p
    a %= p
    if jacobi(a, p) == -1:
        return None
    elif jacobi(a, p) == 0:
        return 0
    r, k = desc(p)
    z = pow(a, r, p)
    x = pow(a, (r + 1) // 2, p)
    if z == 1:
        return x
    s = k - 1
    g = gen_raiz_prim(p)
    while s != 0:
        disc = pow(z, 2 ** (s - 1), p)
        if disc == (p - 1):
            exp = r * (2 ** (k - s - 1))
            t = pow(g, exp, p)
            x *= t
            x %= p
            z *= t ** 2
            z %= p
        s -= 1
    return x

# --------- Resolvemos el caso general

def veces_div_p(a, p):
    k = 0
    while a % p == 0:
        a //= p
        k += 1
    return k

def sqrt_mod(a, p, n):
    x = sqrt_mod1(a, p)
    # caso base
    if n == 1:
        return x
    # a puede no ser coprimo con p: en ese caso, para que haya solución,
    # si a = b * (p ** k) con k >= 1 y p no divide a b,
    # entonces k tiene que ser par y b residuo cuadrático módulo p ** n
    k = veces_div_p(a, p)
    b = a // (p ** k)
    if k != 0:
        if (k % 2) != 0 or jacobi(b, p ** n) != 1:
            return None
        else:
            y = sqrt_mod(b, p, n)
            return (y * (p ** (k // 2))) % (p ** n)
    else:
        # x^2 = a (p ** n) tiene solución si y solo si x^2 = a (p) tiene solución
        if x is None:
            return None
        s = 1
        while s < n:
            t = ((x ** 2) - a) // (p ** s)
            inverso = pow(-2 * x, -1, p)
            l = (t * inverso) % p
            x += l * (p ** s)
            s += 1
        return x

print([sqrt_mod(a, 17, 1) for a in range(17)])
print(sqrt_mod(3, 28091881, 1))
print(sqrt_mod(3, 28091881, 4))
print(sqrt_mod(167042, 17, 7))
print(sqrt_mod(250563, 17, 5))
