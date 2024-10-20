'''
GRUPO: Beñat Pérez de Arenaza Eizaguirre, Mario Calvarro Marines y Diego Ostos Arellano
Comenzamos en la función principal pollard a la que le llega un número n y un valor B que es el límite superior para determinar 
si p-1 es B-suave. El objetivo de esta función es encontrar un primo de n utilizando el método p - 1 de Pollard. 
Comenzamos con la selección de un número aleatorio a entre 2 y n-1, después, calculamos el máximo común divisor(g) entre a y n mediante
la función gcd(a,n), si g es un divisor no trivial entonces hemos terminado. Si no, hay que calcular y = gcd(a^β -1, n) esto nos devuelve y, 
si y no es un divisor trivial de n, entonces tenemos nuestra solución y hemos terminado. No obstante, no calculamos esto directamente, sino que 
para ello usamos la función expo(a, B, n), que genera todos los primos hasta B utilizando el algoritmo de la criba de aristóteles (primos_hasta(B)) 
y con esto, cada primo p, se eleva a la potencia de p(con el objetivo de calcular p^(log p N)). Repitiendo esto para potencias de p hasta que 
excedan a n(es decir elevamos p^(log p N) pero sin usar la función log), realizando las operaciones módulo n. 
'''

import random

def gcd(a, b):
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

def expo(a, B, n):
    """ Calcula a^beta pero haciendo modulo n en cada paso para todos los primos p <= B """
    primos = primos_hasta(B)
    a_beta = a
    for p in primos:
        potencia = p
        a_beta = exp_rapida(a_beta, p, n)
        while potencia < n:
            potencia *= p
            a_beta = exp_rapida(a_beta, p, n)
    return a_beta

def pollard_pm1(n, B):
    """ Implementación del método p-1 de Pollard para factorizar N """
    while True:
        # Paso 1: Elegir a aleatoriamente en [1, N-1]
        a = random.randint(2, n - 1)

        # Paso 2: Calcular gcd(a, N)
        g = gcd(a, n)

        if g > 1 and g < n:
            # Si gcd(a, N) es un divisor no trivial, hemos encontrado un factor
            return g

        # Paso 3: Calcular y = gcd(a^β - 1, N)
        a_beta = expo(a, B, n)
        y = gcd(a_beta - 1, n)

        # Si y es un divisor no trivial de N, lo devolvemos
        if 1 < y < n:
            return y

        # Si y == N, repetimos el proceso
        if y == n:
            continue
