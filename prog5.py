'''
GRUPO: Beñat Pérez de Arenaza Eizaguirre, Mario Calvarro Marines y Diego Ostos Arellano
El objetivo del código es obtener una solución de x^2 ≡ a (mod p^n) donde p es primo impar y la solución que se devuelve debe estar en el intervalo [0,p^n-1].
A lo largo del código usamos varias funciones auxiliares: gcd(a,b) función que recibe dos números y se encarga de devolver el máximo común divisor. jacobi(a,N) es usada para verificar si existe una raíz cuadrada de a mod p. 
Si devuelve -1, no existe raíz cuadrada, si devuelve 1 existe una raíz cuadrada. Si devuelve 0 entonces a es múltiplo de p. La función desc recibe un p y descomponemos p-1 en r * 2^k. 
La función gen_raiz_prim(p) recibe un número p y se encarga de manera aleatoria generar una raíz primitiva módulo p. Otra función auxiliar mas es sqrt_mod1(a,p) que tiene como propósito calcular la raíz cuadrada de a módulo p cuando p es primo,
esto lo hacemos aplicando el algoritmo de Tonelli-Shanks. veces_div_p(a,p) nos devuelve cuantas veces p divide a, es decir la multiplicidad de la raíz p para a. Por último, la función sqrt_mod(a,p,n) que usamos a modo de main y se encarga 
de calcular la solución de x^2 ≡ a (mod p^n). Comenzamos comprobando el símbolo de jacobi para comprobar si existe una raíz cuadrada de a módulo p, si no existe devolvemos None, en caso contrario, continuamos. 
En caso de que n sea 1, ya tenemos la raíz calculada gracias a sqrt_mod1. Mientras que si n > 1, continuamos. Lo que hacemos es ir aplicando múltiples veces el algoritmo concretamente las k-veces que son devueltas por
veces_div_p hasta que k sea igual a 0, hacemos vamos aplicando el algoritmo de una manera recursiva. Hasta llegar al caso base, cuando k = 0, en el recorremos n veces el bucle, lo que hacmoes es calcular la raíz cuadrada de un número a módulo p^n.
Esto es necesario para poder extener la solución desde p hasta p^n. Finalmente, devolvemos el valor de la raíz calculada. 

'''

import random

def gcd(a, b):
    while b > 0:
        a = b
        b = a % b
    return a


def jacobi(a, N):
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
            a = N % a
            N = a
    return m if a != 0 else 0


def desc(p):
    r = p - 1
    k = 0
    while r % 2 == 0:
        k += 1
        r //= 2
    return r, k


def gen_raiz_prim(p):
    while True:
        a = random.randint(2, p - 1)
        if pow(a, (p - 1) // 2, p) != 1:
            return a


def sqrt_mod1(a, p):
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
    k = 0
    while a % p == 0:
        a //= p
        k += 1
    return k


def sqrt_mod(a, p, n):
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
'''
print([sqrt_mod(a, 17, 1) for a in range(17)])
print(sqrt_mod(3, 28091881, 1))
print(sqrt_mod(3, 28091881, 4))
print(sqrt_mod(167042, 17, 7))
print(sqrt_mod(250563, 17, 5))
'''