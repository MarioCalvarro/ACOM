''''
GRUPO: Beñat Pérez de Arenaza Eizaguirre, Mario Calvarro Marines y Diego Ostos Arellano
La función  f(b) = 9b^3 - 225b^2 + 2024b - 5555 - a es estríctamente creciente, puesto que f(b)' = 27b^2 - 450b +2024 > 0 
para todo b en R, ya que no tiene raices reales y f(0)' = 2024 > 0. Además, como el limite de f cuando b tiende a -infinito es -infinito
y el límite de f cuando b tiendo a infinito es infinito, sabemos que para todo valor de a en R exite una única raíz real.


Para el caso en Z, hemos utilizado el método de la bisección. Para ello, debemos encontrar un intervalo [c,d] con c, d en Z
tal que f(c) <= 0 y f(d) >= 0. Para lograr esto, hemos comenzado con c = - |a| y d = |a|. De esta forma, mientras no se cumpla
la condición recién mencionada sumamos una unidad (para evitar problemas en caso de que a = 0) y multiplicamos por dos
c o d en función de cuál de los dos no cumpla la condición de f(c) <= 0 o f(d) >= 0.
(Esto lo hacemos con la función intervalos(a,tipo), donde tipo es un parámetro para diferenciar entre el caso Z y Q).


Una vez encontrado un intervalo tal, procedemos a aplicar el método. Básicamente, calculamos el punto medio i entero de c y d
y evaluamos la función en dicho punto. Como siempre se va a dar que f(c) <= 0 y f(d) >= 0, si f(i) >= 0, sabemos 
que la solución estará en el intervalo [c, i]. En caso de que si f(i) <= 0, la solución estará en el intervalo [i, d].
Este proceso no termina hasta que alguno de los extremos de los intervalos sea una raíz de la función f, o hasta
llegar a un intervalo de longitud 1. Como todos los valores de c y d son enteros, sabemos que estas condiciones de 
parada garantizan encontrar la solución o en su defecto, garantizan que no existe una raíz entera, puesto que
si llegamos a un intervalo de longitud 1 y ninguno de los extremos es raíz de la función, como ambos dos son enteros
no existe ningún entero entre estos últimos, por lo que la raíz no puede ser entera.
(Esto lo hacemos con el método biseccion(c,d,a,tipo)).


Para el caso en Q, el proceso es el mismo hasta llegar a una solución entera o hasta llegar a un intervalo de longitud 1.
En este último caso, todavía no podemos asegurar que no existe raíz racional, pero en caso de existir, debe estar en
[c, c+1]. Por otro lado, por el teorema de la raíz racional, sabemos que si f tiene una raíz racional debe ser de la forma
p/q con q en {1, -1, 3, -3, 9, -9}. No obstante, como habremos encontrado un intervalo de números positivos o negativos 
(es decir, como c es entero, si c < 0 entonces c+1 <= 0, y si c >= 0, entonces c +1 >= 0), reducimos dicho conjunto a {1,3,9}.
Además, como q = 1 se relaciona con los números enteros, ya hemos visto que si hemos llegado hasta este punto sabemos
que f no tiene raíz entera, con lo que nos queda el conjunto {3,9} para q. Esto se traduce en que si la solución debe estar
entre c y c+1, es decir entre 9c/9 y 9(c+1)/9, si existe raíz racional debe ser uno de los siguientes candidatos:
(9c+1)/9, (9c+2)/9, (9c+3)/9, (9c+4)/9, (9c+5)/9, (9c+6)/9, (9c+7)/9, (9c+8)/9. Si nos damos cuenta, el caso de q =3 
también lo estamos incluyendo en (9c+3)/9, (9c+6)/9, por lo que con esos 8 números racionales entre c y c+1 estaríamos
cubriendo todas las posibles raíces racionales de la función f. Por todo ello, evaluaremos la función en estos 8
números racinonales y en caso de que ninguno de ellos sea raíz de f, podemos garantizar que f no tiene raíz racional,
y que por tanto a no es solución racional de a = 9b^3 - 225b^2 + 2024b - 5555. No obstante, para evitar problemas de redondeo, 
hemos decidido hacer lo siguiente. Si tomamos y = 9x, entonces obtenemos la funcion f(x)=f(y/9)= 9^(-2)b^3-225*9^(-2)*b^2+2024*9^(-1)*b-5555-a. Como
estamos buscando la raíz, la siguiente función comparte raíz con la anterior : g(x)= 9^2*f(x) = b^3 - 225b^2 + 2024*9*b - 81(5555+a).
(Para esto último es justamente la función evaluar2(a,b)).
(Todo esto lo comprobamos en la función probar_9)

'''

def evaluar(a, b):
    return ((9*b - 225)*b + 2024)*b - 5555 - a
def evaluar2(a, b):
    return b**3 - 225*b**2 +2024*9*b - 81*(5555+a)

def probar_9(c, a):
    for i in range(1,9):
        if evaluar2(a, 9*c+i) == 0:
            return True
    return False

def biseccion(c, d, a, tipo):
    while True:
        i = (c+d)//2
        evC = evaluar(a, c)
        evD = evaluar(a, d)
        j = evaluar(a, i)
        if j == 0 or evC == 0 or evD == 0:
            return True
        elif evC > 0 or evD < 0 or d-c == 1 and tipo == "z":
            return False
        elif tipo == "q" and d-c == 1:
            return probar_9(c, a)
        elif j < 0:
            c = i
        else:
            d = i

def intervalos(a, tipo):
    c = -abs(a)
    d = abs(a)
    evC = evaluar(a, c)
    evD = evaluar(a, d)
    while evC > 0 or evD < 0:
        if evC > 0:
            c = (c-1)*2
            evC = evaluar(a, c)
        else:
            d = (d+1)*2
            evD = evaluar(a, d)

    return biseccion(c, d, a, tipo)

def tiene_sol_z(a):
    return intervalos(a, "z")

def tiene_sol_q(a):
    return intervalos(a, "q")

