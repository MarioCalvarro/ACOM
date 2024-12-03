# Mario Calvarro, Beñat Pérez, Diego Ostos
##############
# Apartado 1 #
##############
'''
En este apartado, debemos implementar dos funciones mul_pol_mod(f,g,p) y mul_ss_pol_mod(f,g,k,p). La primera función recibe f y g que son polinomios pertenecientes a
(Z/pZ)[x] nos devuelve el polinomio fruto del producto entre f y g. Para la implementación de esta primera función nos hemos apoyado en la segunda función definida 
a continuación. La segunda función realiza el producto de los polinomios f y g en el anillo (Z/pZ)[x](x^(s^k) + 1), devolviendo un resultado respresentado por sus 2^k coeficientes
en orden creciente de grado del monomio. Dividimos los polinomios en representación con u y v mediante la función split. Aplicamos el algoritmo de Cooley-Tuckey para calcular DFT_n, 
transformando así los bloques. Realizamosel producto en el dominio transformado y posteriormente, se usa la transformanda inversa para volver al dominio original. Finalmente, 
unimos los distintos bloques mediante la función join.
'''
def rotar(f, t):
    """
    Realiza la rotación de los elementos del polinomio f en t posiciones, si un valor se sale de rango, entonces se invierte cambiando su rango
    """
    if t == 0:
        return f

    n = len(f)
    res = [0] * n

    for i in range(n):
        if i + t >= n or i + t < 0:
            res[( i + t ) % n] = -f[i]
        else:
            res[( i + t ) % n] = f[i]

    return res

def sumar_vectores(f, g, p):
    """
    Dados dos vectores f y g realizamos la suma de ambos en módulo p
    """
    h = []
    for i in range(len(f)):
        h += [(f[i] + g[i]) % p]
    return h

def restar_vectores(f, g, p):
    """
    Dados dos vectores f y g realizamos la resta de ambos en módulo p
    """
    h = []
    for i in range(len(f)):
        h += [(f[i] - g[i]) % p]
    return h
 
def fft_adaptado(f, xi, p):
    """
    Algoritmo de Cooley-Tuckey para calcular DFT_n(p) donde n = len(p) debe ser una potencia de 2 y xi debe ser una raíz n-ésima primitiva de la unidad
    """
    n = len(f)
    if n == 1:
        return f

    p_even = [0] * (n//2)
    p_odd  = [0] * (n//2)
    for i in range(n//2):
        p_even[i] = f[2*i]
        p_odd[i]  = f[2*i+1]

    a_even = fft_adaptado(p_even, xi*2, p)
    a_odd  = fft_adaptado(p_odd,  xi*2, p)

    a = [0] * n
    for i in range(n//2):
        rotar_impar = rotar(a_odd[i], i*xi)

        a[i] = sumar_vectores(a_even[i], rotar_impar, p)
        a[i + n // 2] = restar_vectores(a_even[i], rotar_impar, p)

    return a

def ifft_adaptado(a, xi, p):
    """
    Esta función calcula la transformada inversa utilizando que D(xi)^(-1) = 1/n * D(1/xi)
    """
    n1 = len(a[0])
    n2 = len(a)
    transformado = fft_adaptado(a, -xi, p)
    inversa = pow(n2, -1, p)

    for i in range(n2):
        for j in range(n1):
            transformado[i][j] = (transformado[i][j] * inversa) % p

    return transformado

def split(f, k):
    """
    Divide un polinomio en representación con u y v.
    """
    k2 = (k + 1) // 2
    k1 = k - k2

    n1 = 2**k1
    n2 = 2**k2
    res = []
    for i in range(0, n2):
        res += [f[i*n1 : i*n1 + n1] + [0]*n1]

    return res

def join(f_split, p):
    """
    Reconstruye un polinomio desde su representación con u y v.
    """
    n2 = len(f_split)
    n1 = len(f_split[0])//2

    pol = restar_vectores(f_split[0][0:n1], f_split[n2-1][n1:(2*n1)], p)
    for i in range(0,n2-1):
        pol   = pol + sumar_vectores(f_split[i][n1:(2*n1)], f_split[i+1][0:n1], p)

    return pol

def negaproducto(f, n1, n2):
    """
    Realiza el producto w·f de la negaconvolucion. En este caso, son
    simplemente rotaciones de los elementos de f
    """
    t = (2*n1) // n2
    res = [f[0]]

    for i in range (1, n2):
        res += [rotar(f[i], i*t)]

    return res

def inv_negaproducto(f, n1, n2):
    """
    Realiza el producto w·f de la negaconvolucion. En este caso, son
    simplemente rotaciones inversas de los elementos de f
    """
    t = (2*n1) // n2
    res = [f[0]]

    for i in range (1, n2):
        res += [rotar(f[i], -(i*t))]

    return res

def negaconvolucion(f, g, k1, p):
    """
    Recibe dos polinomios f y g, un entero k1 y el módulo del espacio en el que estamos. La usamos para realizar la multiplicación eficiente de los polinomios diviendolos en bloques. 
    """
    n1 = len(f[0])//2
    n2 = len(f)

    term1 = fft_adaptado(negaproducto(f, n1, n2), (4*n1)//n2, p)
    term2 = fft_adaptado(negaproducto(g, n1, n2), (4*n1)//n2, p)

    prod = []
    for i in range(n2):
        prod += [mul_ss_pol_mod(term1[i], term2[i], k1 + 1, p)]

    transformado_inv = ifft_adaptado(prod, (4*n1)//n2, p)

    return inv_negaproducto(transformado_inv, n1, n2)

def mul_ss_pol_mod(f, g, k, p):
    """
    Recibe dos polinomios f y g, k para los coeficientes que se usan para la representación y p define el módulo en el que estamos.
    Devolvemos el resultado del producto de ambos polinomios representado por la lista de sus 2^k coeficientes.
    """
    if f == [] or g == []:
        return []

    # Casos base:
    if k == 0:
        return [(f[0]*g[0]) % p]

    if k == 1:
        return [(f[0]*g[0] - f[1]*g[1])%p, (f[0]*g[1] + f[1]*g[0]) % p]

    if k == 2:
        a_0 = (f[0]*g[0] - f[1]*g[3] - f[2]*g[2] - f[3]*g[1]) % p
        a_1 = (f[0]*g[1] + f[1]*g[0] - f[2]*g[3] - f[3]*g[2]) % p
        a_2 = (f[0]*g[2] + f[1]*g[1] + f[2]*g[0] - f[3]*g[3]) % p
        a_3 = (f[0]*g[3] + f[1]*g[2] + f[2]*g[1] + f[3]*g[0]) % p
        return [a_0, a_1, a_2, a_3]

    # Caso recursivo:
    k2 = (k + 1) // 2
    k1 = k - k2

    f_class = split(f, k)
    g_class = split(g, k)

    h_class = negaconvolucion(f_class, g_class, k1, p)

    return join(h_class, p)

def eliminar_ceros(v):
    """
    Recibe un polinomio de grado n y elimina los ceros desde la posición n - 1 hasta el primer número no nulo.
    """
    n = len(v)
    while n >= 1 and v[n-1] == 0:
        n -= 1
    del v[n:]
    return v

def mul_pol_mod(f, g, p):
    """
    Recibe dos polinomios f y g en módulo p
    Devuelve como resultado el polinomio producto de la multiplicación de ambos polinomios aplicando la función mul_ss_pol_mod en k, obtenida por medio 
    de dividir m entre dos hasta que sea cero, "cuantas veces" sea m divisible por 2. 
    """
    if f == [] or g == []:
        return []

    m = len(f) + len(g)
    k = 0
    while m != 0:
        m //= 2
        k += 1

    # Añadimos los ceros necesarios
    f = f + [0]*(2**k - len(f))
    g = g + [0]*(2**k - len(g))

    h = eliminar_ceros(mul_ss_pol_mod(f, g, k, p))
    return h

##############
# Apartado 2 #
##############
"""
El objetivo de este apartado es realizar la implementación de tres funciones (explicadas posteriormente)
"""

def mul_lo_toep_mod(v, a, p, n):
    """
    Calcula el producto de L(v) * a, donde L(v) es una matriz de Toeplitz inferior cuya primera columna es v y a es un vector. 
    Esto se obtiene multiplicando el polinomio v por a y rellenando con ceros hasta obtener un tamaño n.
    """
    aux = mul_pol_mod(v, a, p)

    if len(aux) < n:
        aux += [0] * (n - len(aux))

    return aux[:n]

def mul_hi_toep_mod(v, a, p, n):
    """
    Calcula el producto de U(v) * a, donde U(v) es una matriz de Toeplitz superior cuya primera fila es v y el vector a. 
    Para realizar esto invertimos el orden de a, realizamos la multiplicación con el polinomio v y finalmente, tomamos los n primeros coeficientes. 
    En caso de que nos quede un polinomio de grado menor que n, rellenamos con ceros
    """
    a_aux = a[::-1]
    aux = mul_pol_mod(v, a_aux, p)
    aux_2 = aux[:n]

    if len(aux_2) < n:
        aux_2 += [0] * (n - len(aux_2))

    return aux_2[::-1]

def mul_toep_mod(v, w, a, p, n):
    """
    Calcula el producto de T(v,w) * a, donde T(v,w) es una matriz de Toeplitz general cuya primera columna es v y su primera filaes w
    Para realizar esto, calculamos U(v) * a y L(w) * a, rellenemos con ceros hasta que ambas tenga mismo grado n. Y terminamos sumando ambos términos.
    """
    term1 = mul_hi_toep_mod(v, a, p, n)
    term2 = [0] + (mul_pol_mod(w, a[0:n-1], p))[0:n-1]
    
    if len(term1) < n:
        term1 += [0] * (n - len(term1))
    if len(term2) < n:
        term2 += [0] * (n - len(term2))

    return sumar_vectores(term1, term2, p)

##############
# Apartado 3 #
##############
"""
El objetivo de este apartado es calcular dos funciones, comentadas a continuación
"""
def mul_inv_hi_toep_mod(v, a, p, n):
    """
    Calcula el producto U^-1(v) * a, donde U(v) es una matriz de Toeplitz superior cuya cuya primera fila es v y el vector a. 
    Casos base: 
        Si n = 1: devolvemos producto del inverso v[0] y a[0]
        Si n = 2: aplicamos el inverso a ambos términos del resultado. 
    Caso recusivo: 
        Dividimos v y a en dos mitades, resolvemos de manera recursiva los productos con las submatrices superiores de menor tamaño. Nos apoyamos en la función mul_toep_mod para obtener los resultados.
        Finalmente, unimos los resultados para obtener los resultados.
    Si n no es par, extendemos v y a con ceros.
    """
    if n == 1:
        inverso = pow(v[0], -1, p)
        return [(inverso * a[0]) % p]
    elif n == 2:
        inverso = pow(v[0], -1, p)
        return [(inverso * a[0] - v[1] * pow(inverso, 2, p) * a[1]) % p,
                (inverso * a[1]) % p]

    if n % 2 == 0:
        t1 = v[n//2:n]
        t2 = v[n//2-1:0:-1]
        mul_uinv_a1 = mul_inv_hi_toep_mod(v[0:n//2], a[0:n//2], p, n//2)
        mul_uinv_a2 = mul_inv_hi_toep_mod(v[0:n//2], a[n//2:n], p, n//2)

        mul_t = mul_toep_mod(t1, t2, mul_uinv_a2, p, n//2)

        mul_uinv_t_a2 = mul_inv_hi_toep_mod(v[0:n//2], mul_t, p, n//2)
        return restar_vectores(mul_uinv_a1, mul_uinv_t_a2, p) + mul_uinv_a2
    else:
        return (mul_inv_hi_toep_mod(v + [0], a + [0], p, n + 1))[0:n]

def mul_inv_lo_toep_mod(v, a, p, n):
    """
    Calculamos L^-1(v) * a, donde L(v) es la matriz de Toeplitz inferiro cuya primera columna es v y el vector a. Usamos el complemento de Schur para ver que la fórmula es simétrica
    Casos base: 
        Si n = 1: devolvemos producto del inverso v[0] y a[0]
        Si n = 2: aplicamos el inverso a ambos términos del resultado.
    Caso recusrivo: 
        Dividimos v y a en dos submatrices, así calculamos recursivamente las submatrices inferiores de menor tamaño. Nos apoyamos en la función mul_toep_mod para obtener los resultados.
        Finalmente, unimos los resultados para obtener los resultados.
    Si n no es par, extendemos v y a con ceros. 
    """
    if n == 1:
        inverso = pow(v[0], -1, p)
        return [(inverso * a[0]) % p]
    elif n == 2:
        inverso = pow(v[0], -1, p)
        return [(inverso * a[0]) % p,
                (inverso * a[1] - v[1] * pow(inverso, 2, p) * a[0]) % p]

    if n % 2 == 0:
        t1 = v[n//2:0:-1]
        t2 = v[n//2+1:n]
        mul_linv_a1 = mul_inv_lo_toep_mod(v[0:n//2], a[0:n//2], p, n//2)
        mul_linv_a2 = mul_inv_lo_toep_mod(v[0:n//2], a[n//2:n], p, n//2)

        mul_t = mul_toep_mod(t1, t2, mul_linv_a1, p, n//2)

        mul_linv_t_a1 = mul_inv_lo_toep_mod(v[0:n//2], mul_t, p, n//2)
        return mul_linv_a1 + restar_vectores(mul_linv_a2, mul_linv_t_a1, p)
    else:
        return (mul_inv_lo_toep_mod(v + [0], a + [0], p, n + 1))[0:n]

##############
# Apartado 4 #
##############
"""
El objetivo de esta parte es dados dos polinomios f y g pertenenecientes a Z/pZ de grado n y m tal que n >= m, obtenemos la división con resto, donde el cociente es
de grado n - m = deg(q) y deg(r) < m. En resumen, f = q * g + r.
La función divmod_pol_mod(f,g,p) que devuelve q y r
"""
def divmod_pol_mod(f, g, p):
    """
    Calcula la división de dos polinomios f y g en un cuerpo Z/pZ, devolviendo el cociente y el resto de la división bajo módulo p. 
    Si el grado de f es menor que el grado de g, devolvemos [] y f.
    Si el grado de f es mayor o igual que el grado de g, utilizamos las matrices de Toeplitz superiores para resolver el sistema lineal asociado a la división y calculando así el cociente mediante la función
    mul_inv_hi_toep_mod. Para calcular el resto, aplicamos la función restar vectores con f y el resultado de multiplicar el cociente por el divisor. 
    En ambos calculos, finalizamos con la eliminación de los ceros innecesarios.
    """
    n = len(f) - 1
    m = len(g) - 1
    # Caso deg(f) < deg(g)
    if n < m:
        return [], f

    # Caso deg(f) >= deg(g)
    final_g = max(-1, 2*m - n - 2)
    matriz_u = []
    if final_g == -1:
        matriz_u = g[m::-1]
        matriz_u += [0] * (n - m + 1 - len(matriz_u))
    else:
        matriz_u = g[m:final_g:-1]

    vect_f = f[m:n+1]
    
    vect_q = eliminar_ceros(mul_inv_hi_toep_mod(matriz_u, vect_f, p, n - m + 1))

    vect_r = eliminar_ceros(restar_vectores(f, mul_pol_mod(vect_q, g, p), p))

    return vect_q, vect_r

