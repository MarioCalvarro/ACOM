from prog6 import *
#Casos de prueba para la Parte 4
print("Caso de prueba 1")
f = [5, 3, 4]
g = [3]
p = 5
print(divmod_pol_mod(f, g, p) == ([0, 1, 3], []))
print("Caso de prueba 2")
f = [1, 2, 4, 3]
g = [1, 2]
p = 7
print(divmod_pol_mod(f, g, p) == ([3, 3, 5], [5]))
print("Caso de prueba 3")
f = [1, 4, 3, 2]
g = [1, 1, 1]
p = 5
print(divmod_pol_mod(f, g, p) == ([1,2],[0,1]))
print("Caso de prueba 4")
f = [1, 2, 3, 4]
g = [1, 1]
p = 7
print(divmod_pol_mod(f, g, p) == ([3,6,4],[5]))
print("Caso de prueba 5")
f = [1, 1, 1, 1,1]
g = [1, 1, 1,1]
p = 7
print(divmod_pol_mod(f, g, p) == ([0, 1],[1]))
print("Caso de prueba 6")
f = [1, 4, 6]
g = [1, 2]
p = 5
print(divmod_pol_mod(f, g, p) == ([3, 3],[3]))
print("Caso de prueba 7")
f = [3, 1]
g = [1, 2, 4]
p = 7
print(divmod_pol_mod(f, g, p) ==([],[3, 1]))
