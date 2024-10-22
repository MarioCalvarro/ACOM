def tonelli_shanks(a, p):
    """Encuentra x tal que x^2 ≡ a (mod p) usando Tonelli-Shanks.
    Devuelve None si no hay solución."""
    if pow(a, (p - 1) // 2, p) != 1:
        return None  # No hay solución

    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)

    # Tonelli-Shanks general
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1

    z = 2
    while pow(z, (p - 1) // 2, p) == 1:
        z += 1

    m = s
    c = pow(z, q, p)
    t = pow(a, q, p)
    r = pow(a, (q + 1) // 2, p)

    while t != 1:
        i = 0
        temp = t
        while temp != 1:
            temp = pow(temp, 2, p)
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = pow(b, 2, p)
        t = (t * c) % p
        r = (r * b) % p

    return r

def hensel_lift(a, p, n, x0):
    """Aplica el método de Hensel para levantar la solución de x^2 ≡ a (mod p^n) a (mod p^(n+1))"""
    x = x0
    pn = p
    for i in range(1, n):
        if (x * x - a) % pn != 0:
            return None  # No hay solución
        f_prime = (2 * x) % p
        if f_prime == 0:
            return None  # No hay solución
        x = (x + ((a - x * x) // pn) * pow(f_prime, -1, p)) % (pn * p)
        pn *= p
    return x

def sqrt_mod(a, p, n):
    """Calcula una solución para x^2 ≡ a (mod p^n), donde p es un primo impar y n ≥ 1."""
    if n == 1:
        return tonelli_shanks(a, p)

    # Primero encuentra la solución módulo p
    x0 = tonelli_shanks(a, p)
    if x0 is None:
        return None

    # Aplica el método de Hensel para elevar la solución a módulo p^n
    return hensel_lift(a, p, n, x0)
