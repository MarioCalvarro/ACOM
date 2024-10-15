import random
import math

def gcd(x, y):
    """ Calcula el máximo común divisor (gcd) de x e y """
    x = abs (x)
    y = abs (y)
    xespar = x %2 == 0
    yespar = y %2 == 0
    if x == 0: # caso base : gcd (0 , y )= y
        m = y
    elif y == 0: # caso base : gcd (x ,0)= x
        m = x
    elif xespar and yespar :
        m = 2 * gcd( x //2 , y //2)
    elif xespar :
        m = gcd( x //2 , y )
    elif yespar :
        m = gcd(x , y //2)
    elif x > y :
        m = gcd(y , (x - y )//2)
    else :
        m = gcd(x , (y - x )//2)
    return m

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
    primes = primos_hasta(B)
    beta_val = a
    for p in primes:
        potencia = p
        beta_val = (beta_val ^ p) % n
        while potencia < n:
            potencia *= p
            beta_val = (beta_val ^ p) % n
    return beta_val

def pollard_pm1(n, B):
    """ Implementación del método p-1 de Pollard para factorizar N """
    # Paso 1: Calcular β
    ##beta_val = beta(B, n)

    while True:
        # Paso 2: Elegir a aleatoriamente en [1, N-1]
        a = random.randint(2, n - 1)

        # Paso 3: Calcular gcd(a, N)
        g = gcd(a, n)

        if g > 1 and g < n:
            # Si gcd(a, N) es un divisor no trivial, hemos encontrado un factor
            return g

        # Paso 4: Calcular y = gcd(a^β - 1, N)
        a_beta = expo(a, B, n)
        ##a_beta = exp_rapida(a, beta_val, n)
        y = gcd(a_beta - 1, n)

        # Si y es un divisor no trivial de N, lo devolvemos
        if 1 < y < n:
            return y

        # Si y == N, repetimos el proceso
        if y == n:
            continue

#pollard_pm1(266941026756980432674273004285227481087654764302423743652775141424617302734658657044184670961906473380670375458124002139355389903179856802994870173312715601752745436095944859694011281380544008461966716092656817602170627934379,1000000)
print(pollard_pm1(6762691917284635893988303,100)) # 263381007359
print(pollard_pm1(119495316416051485815407352901774402921283,1000)) # 45193190000498511563
print(pollard_pm1(58771624767109592190294695436237320281954733562633238951,10000)) # 455299501714901028162586319
print(pollard_pm1(266941026756980432674273004285227481087654764302423743652775141424617302734658657044184670961906473380670375458124002139355389903179856802994870173312715601752745436095944859694011281380544008461966716092656817602170627934379,1000000))# 369533287075647816702693047063862986159754020633236293082394185958849017640739268469321177775535208522441080367