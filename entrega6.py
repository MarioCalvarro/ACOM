# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:41:47 2023

@author: Juan Vallaure Lozano y Carlos Antonio Sánchez-Beato Alonso
"""

# Entrega 6 (problema 7.7):

#IMPLEMENTACIÓN DE mult_ss_mod
# Dados f y g polinomios en (Z/pZ)[X]/(X^n + 1) con gr(f)=gr(g)=n-1 y n=2^k.
# Visualizamos estos polinomios en A = {(Z/pZ)[u]/(u^2n1 + 1)}[z] / (z^n2 + 1),
# con los n1 y n2 adecuados tales que n = n1n2 y 2n1>=n2.

# En este anillo se cumplen las condiciones para realizar una negaconvolución.
# Sabemos que aquí [u^(2n1/n2)] es una raiz 2n2-ésima.
# Obtenemos inmediatamente que [u^4(n1/n2)] es una raiz n2-ésima.
# Multiplicar por estos elementos en A será rotar los coeficientes (polinomios).


# Función auxiliar 1: rota los elementos de una lista t casillas
def shift (f, t):
    if t == 0:
        return f

    l = len(f)
    aux = [0]*l

    for j in range(l):
        # la rotación tiene en cuenta que [u^n2] = [-1]
        if j+t >= l or j+t<0:
            aux[ ( j + t ) % l] = -f[j]

        else:
            aux[ ( j + t ) % l] = f[j]

    return aux


# Función auxiliar 2: hace la operación w·f de la negaconvolución (2n1(n2-1) ops.)
def nega_shift (f, n1, n2):
    t = (2*n1)//n2
    nega = [f[0]]

    for i in range (1, n2):
        # w = (1, [u^(2n1/n2)], [u^(2n1/n2)]^2, ..., [u^(2n1/n2)^(n2-1)])
        aux = shift(f[i], i*t)
        nega += [aux]

    return nega

# Función auxiliar 3: hace la operación w^-1·f de la negaconvolución (2n1(n2-1) ops.).
# Obsérvese que w^-1 corresponde a ejecutar las rotaciones opuestas.
def anti_nega_shift (f, n1, n2):
    t = (2*n1)//n2
    nega = [f[0]]

    for i in range (1, n2):
        aux = shift(f[i],-(i*t))
        nega += [aux]

    return nega

# Función auxiliar 4: suma modular (p) de vectores de misma longitud
def suma_vectores_mod(a,b,p):
    l = []
    for i in range(len(a)):
        l += [(a[i]+b[i])%p]
    return l

# Función auxiliar 5: resta modular (p) de vectores de misma longitud
def resta_vectores_mod(a,b,p):
    l = []
    for i in range(len(a)):
        l += [(a[i]-b[i])%p]
    return l

# Función auxiliar 6: transformada de fourier adaptada a los elementos de A.
# Hemos cambiado xi por j, que corresponde al salto de la rotación asociada.
def fft (f, j, p):
    n2 = len(f)
    if n2 == 1:
        return f

    f_even = [0]*(n2//2)
    f_odd  = [0]*(n2//2)

    for i in range (n2//2):
        f_even[i] = f[2*i]
        f_odd[i]  = f[2*i + 1]

    # xi**2 es como rotar dos veces, ie, 2*j
    a_even = fft(f_even, 2*j, p)
    a_odd  = fft(f_odd , 2*j, p)
    a = [0] * n2

    for i in range(n2//2):
        # de nuevo, xi**i es como rotar i veces, ie, i*j
        shift_a_odd = shift (a_odd[i], i*j)

        a[i]         = suma_vectores_mod (a_even[i], shift_a_odd, p)
        a[i + n2//2] = resta_vectores_mod(a_even[i], shift_a_odd, p)

    return a


# Función auxiliar 7: inversa de fft también adaptada a los elementos de A.
def ifft(f, j, p):
    n1 = len(f[0])
    n2 = len(f)
    transformada = fft(f, -j, p)  # tomamos -j para rotar en sentido opuesto.
    inverse = pow(n2,-1,p)

    for i in range(n2):
        for j in range(n1):
            transformada[i][j] = (transformada[i][j] * inverse)%p

    return transformada

# Función auxiliar 8: definimos una función que devuelva el resultado a (Z/pZ)[X]/(X^n + 1)
def volver_al_mundo_normal (l, p):
    n2 = len(l)
    n1 = len(l[0])//2
    # En vez de operar exhaustivamente (O(n(n2-1)2n1)), operamos en bloques
    # Los últimos n1 términos de la última lista se desbordan:
    pol = resta_vectores_mod(l[0][0:n1], l[n2-1][n1:(2*n1)], p)
    for i in range(0,n2-1):
        pol   = pol + suma_vectores_mod(l[i][n1:(2*n1)], l[i+1][0:n1], p)
    return pol


# Ya podemos definir la función recursiva
def mult_ss_mod(f, g, k, p):
    if f == [] or g == []:
        return []

    # Casos base:
    if k == 0:
        return [(f[0]*g[0])%p]

    if k == 1:
        return [(f[0]*g[0] - f[1]*g[1])%p, (f[0]*g[1] + f[1]*g[0])%p]

    if k == 2:
        a_0 = (f[0]*g[0] - f[1]*g[3] - f[2]*g[2] - f[3]*g[1])%p
        a_1 = (f[0]*g[1] + f[1]*g[0] - f[2]*g[3] - f[3]*g[2])%p
        a_2 = (f[0]*g[2] + f[1]*g[1] + f[2]*g[0] - f[3]*g[3])%p
        a_3 = (f[0]*g[3] + f[1]*g[2] + f[2]*g[1] + f[3]*g[0])%p
        return [a_0, a_1, a_2, a_3]

    # Distinguimos los k1 (n1) y k2 (n2) de los que hemos ido hablando
    if k%2 == 0:
        k1, k2 = k//2, k//2
    else:
        k1, k2 = (k-1)//2, (k+1)//2

    n1, n2 = 2**k1, 2**k2

    # Transformamos f y g a elementos de A
    f_class = []
    g_class = []

    for i in range(0, n2):
        f_class += [f[i*n1 : i*n1 + n1] + [0]*n1]
        g_class += [g[i*n1 : i*n1 + n1] + [0]*n1]


    # Obtenemos w·f y w·g
    w_f = nega_shift(f_class, n1, n2)
    w_g = nega_shift(g_class, n1, n2)

    # Realizamos DFTn(w·f) y DFTn(w·g)
    w_f_fft = fft(w_f, (4*n1)//n2, p)
    w_g_fft = fft(w_g, (4*n1)//n2 ,p)

    # Multiplicamos recursivamente y coordenada a coordenada DFTn(w·f) y DFTn(w·g).
    producto =  []
    for i in range(n2):
        producto += [mult_ss_mod(w_f_fft[i], w_g_fft[i], k1 + 1, p)]

    # al resultado le aplicamos la inversa de la transformada de fourier
    ifft_producto = ifft(producto, (4*n1)//n2, p)
    # y luego le multiplicamos por w^-1.
    solucion = volver_al_mundo_normal(anti_nega_shift(ifft_producto, n1, n2), p)

    # Finalmente, mandamos el polinomio obtenido de nuevo a (Z/pZ)[X]/(X^n + 1)

    return solucion


# Implementemos ahora mult_pol_mod.
# Bastará encontrar la potencia de 2 que supere a gr(f) + gr(f)
# Pues ese es el grado del polinomio producto f·g
# Luego, ejecutamos  mult_ss_mod(f, g, k, p), usando g y f vistos en (Z/pZ)[X]/(X^n + 1)
# Para eso bastará rellenar con los ceros necesarios
# Finalmente, habrá que quitar los ceros que aparezcan al final del resultado.
def remover_ceros (l):
    n = len(l)
    while n >= 1 and l[n-1] == 0:
        n -= 1
    del l[n:]
    return l

def mult_pol_mod(f, g,p):
    if f == [] or g == []:
        return []

    m = len(f) + len(g)
    k = 0
    while m != 0:
        m //= 2
        k += 1
    f = f + [0]*(2**k - len(f))
    g = g + [0]*(2**k - len(g))
    h = remover_ceros(mult_ss_mod(f, g, k, p))
    return h

primero = [0, 1, 2, 3, 4, 5, 6, 7]
segundo = [0, 10, 20, 30, 40, 50, 60, 70]
k = 3
p = 1009

mult_ss_mod(primero, segundo, k, p)
