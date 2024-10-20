'''
GRUPO: Beñat Pérez de Arenaza Eizaguirre, Mario Calvarro Marines y Diego Ostos Arellano
Comenzamos en la función principal pollard a la que le llega un número n y un valor B que es el límite superior para determinar 
si p-1 es B-suave. El objetivo de esta función es encontrar un primo de n utilizando el método p - 1 de Pollard. 
Comenzamos con la selección de un número aleatorio a entre 2 y n-1, después, calculamos el máximo común divisor(g) entre a y n mediante
la función gcd(a,n), si g es un divisor no trivial entonces hemos terminado. Si no, hay que calcular y = gcd(a^β -1, n) esto nos devuelve y, 
si y no es un divisor trivial de n, entonces tenemos nuestra solución y hemos terminado. No obstante, no calculamos esto dorectamente, sino que para ello usamos la función expo(a, B, n), que genera todos los primos hasta B
utilizando el algoritmo de la criba de aristóteles (primos_hasta(B)) y con esto, cada primo p, se eleva a a la potencia de p(con el objetivo de calcular p^(log p N). Repitiendo esto para
potencias de p hasta que excedan a n(es decir elevamos p^(log p N) pero sin usar la función log), realizando las operaciones módulo n. 
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

# print(pollard_pm1(6762691917284635893988303,100)) # 263381007359
# print(pollard_pm1(119495316416051485815407352901774402921283,1000)) # 45193190000498511563
# print(pollard_pm1(58771624767109592190294695436237320281954733562633238951,10000)) # 455299501714901028162586319
# print(pollard_pm1(266941026756980432674273004285227481087654764302423743652775141424617302734658657044184670961906473380670375458124002139355389903179856802994870173312715601752745436095944859694011281380544008461966716092656817602170627934379,1000000))# 369533287075647816702693047063862986159754020633236293082394185958849017640739268469321177775535208522441080367
print(pollard_pm1(8015110843689733987734493314126305093304900427644362111427136747967674950102502227401067508400666253684402757742551821894809831742615784192918469456560690085286199910126971564784243019741380555941733042571384312328652773515182854067261862538680213978368174296772106253247809802194498297934955888000983604718178338431159811628901253189704085615502336825120630199737136386667975473370270254901099100474927139919509948252709580689468366977236252238872840439164497603293917510631536134556475105354896902911935764652153360729834323397346255158826324977097812765382378557607077261662089097246276601527470070388078314698083342531125763586195801205695074271301766840569856852751612504601753988984656300401973433861641028853770900162658263178765900267910765383841534079073198847564915836063133010526813527,10**7))
