import random
import math

def gcd(a, b):
    """ Calcula el máximo común divisor (gcd) de a y b """
    while b != 0:
        a, b = b, a % b
    return a

def exp_rapida(base, exp, mod):
    result = 1
    base = base % mod
    while exp > 0:
        if exp % 2 == 1:
            result = (result * base) % mod
        exp = exp // 2
        base = (base * base) % mod
    return result

def primos_hasta(n):
    """ Genera los números primos hasta n usando el algoritmo de la criba de Eratóstenes """
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False

    i = 2
    while i * i <= n:
        if sieve[i]:
            for j in range(i*i, n + 1, i):
                sieve[j] = False
        i += 1

    return [x for x in range(2, n + 1) if sieve[x]]

def beta(B, n):
    """ Calcula β como producto de p^ceil(logp(N)) para todos los primos p <= B """
    primes = primos_hasta(B)
    beta_val = 1
    for p in primes:
        beta_val *= exp_rapida(p, int(math.ceil(math.log(B, p))), n)
    return beta_val

def pollard_pm1(n, B):
    """ Implementación del método p-1 de Pollard para factorizar N """
    # Paso 1: Calcular β
    beta_val = beta(B, n)

    while True:
        # Paso 2: Elegir a aleatoriamente en [1, N-1]
        a = random.randint(2, n - 1)

        # Paso 3: Calcular gcd(a, N)
        g = gcd(a, n)

        if g > 1 and g < n:
            # Si gcd(a, N) es un divisor no trivial, hemos encontrado un factor
            return g

        # Paso 4: Calcular y = gcd(a^β - 1, N)
        a_beta = exp_rapida(a, beta_val, n)
        y = gcd(a_beta - 1, n)

        # Si y es un divisor no trivial de N, lo devolvemos
        if 1 < y < n:
            return y

        # Si y == N, repetimos el proceso
        if y == n:
            continue

pollard_pm1(6762691917284635893988303,100)
