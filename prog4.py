import random

def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

def fast_exp(base, exp, mod):
    result = 1
    base = base % mod
    while exp > 0:
        if exp % 2 == 1:
            result = (result * base) % mod
        exp = exp // 2
        base = (base * base) % mod
    return result

def pollard_pm1(n, B):
    a = random.randint(2, n - 2)  # Escoge un valor aleatorio para a
    d = gcd(a, n)  # Comprueba si ya tenemos un factor común

    if d > 1:
        return d  # Si a y n ya comparten un divisor, lo devolvemos

    # Calcula a^M mod n para M = producto de todos los primos ≤ B
    for j in range(2, B + 1):
        a = fast_exp(a, j, n)  # Usamos exponenciación rápida
        d = gcd(a - 1, n)  # Comprueba si a-1 tiene un divisor común con n

        if 1 < d < n:
            return d  # Devuelve el divisor primo encontrado

    return None  # No se encontró factor en el rango dado

# Ejemplo de uso:
n = 8051  # Producto de dos primos p y q
