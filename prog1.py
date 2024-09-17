# Mario Calvarro Marines. Álgebra computacional
# Polinomio: a = 9b^3 - 225b^2 + 2024b - 5555

div_a_n = [1, 3, 9]

def divisores(num):
    if num < 0:
        num = -num

    divs = [i for i in range(1, num + 1) if num % i == 0]
    return divs


def tiene_sol_z(a):
    # El polinomio es estrictamente creciente por lo que solo tendrá una
    # solución real. Se trata de ver si esta es entera.
    term_ind = - (a + 5555)
    divs = divisores(term_ind)

    for b in divs:
        if (b * (b * (9*b - 225) + 2024) + term_ind) == 0:
            return True

    return False

print(tiene_sol_z(5165149559269435554097224863147517932681599350930788324577847807772824643103215545497557949687873329190211893100996934536646532754961536003695719039165489216451850959451311425973420804580551467074979145992385679697710273119451787452893650246060761748772289677435795344815357340032695573403958556958019507982537999231207872612403850126388105936850777470328598547306476039601273718224587269649811422223493357343013581636569090664644393634872128587233785652251717227812667974937))

def tiene_sol_q(a):
    # El polinomio es estrictamente creciente por lo que solo tendrá una
    # solución real. Se trata de ver si esta es racional.
    return False
