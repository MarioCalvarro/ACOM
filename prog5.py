import random

def gcd(a, b):
    """Calcula el máximo común divisor de a y b."""
    while b:
        a, b = b, a % b
    return a


def jacobi(a, N):
    """Calcula el símbolo de Jacobi (a|N)."""
    m = 1
    while a != 0 and a != 1:
        if a < 0:
            if N % 4 == 1:
                m *= -1
            a = -a
        elif a >= N:
            a = a % N
        elif a % 2 == 0:
            if N % 8 in [3, 5]:
                m *= -1
            a //= 2
        else:
            if a % 4 == 3 and N % 4 == 3:
                m *= -1
            a, N = N % a, a
    return m if a != 0 else 0


def desc(p):
    """Descompone p-1 como r * 2^k."""
    r = p - 1
    k = 0
    while r % 2 == 0:
        k += 1
        r //= 2
    return r, k


def gen_raiz_prim(p):
    """Genera una raíz primitiva módulo p."""
    while True:
        a = random.randint(2, p - 1)
        if pow(a, (p - 1) // 2, p) != 1:
            return a


def sqrt_mod1(a, p):
    """Calcula una raíz cuadrada de a módulo primo p usando el algoritmo de Tonelli-Shanks."""
    if jacobi(a, p) == -1:
        return None
    if jacobi(a, p) == 0:
        return 0
    r, k = desc(p)
    z = pow(a, r, p)
    x = pow(a, (r + 1) // 2, p)
    if z == 1:
        return x
    g = gen_raiz_prim(p)
    s = k - 1
    while s != 0:
        disc = pow(z, 2 ** (s - 1), p)
        if disc == p - 1:
            t = pow(g, r * 2 ** (k - s - 1), p)
            x = (x * t) % p
            z = (z * t * t) % p
        s -= 1
    return x


def veces_div_p(a, p):
    """Cuenta cuántas veces p divide a."""
    k = 0
    while a % p == 0:
        a //= p
        k += 1
    return k


def sqrt_mod(a, p, n):
    """Calcula una solución de x^2 ≡ a mod p^n."""
    x = sqrt_mod1(a, p)
    if x is None:
        return None
    if n == 1:
        return x
    k = veces_div_p(a, p)
    b = a // (p ** k)
    if k != 0:
        if k % 2 != 0 or jacobi(b, p ** n) != 1:
            return None
        y = sqrt_mod(b, p, n)
        return (y * (p ** (k // 2))) % (p ** n)
    else:
        s = 1
        while s < n:
            t = ((x ** 2) - a) // (p ** s)
            inv = pow(-2 * x, -1, p)
            l = (t * inv) % p
            x = (x + l * (p ** s)) % (p ** (s + 1))
            s += 1
        return x % (p ** n)

