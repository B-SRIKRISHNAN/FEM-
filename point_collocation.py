import sympy as sym
import math as mt
import numpy as np
import scipy as sy
from matplotlib import pyplot as plt


def residual_calc(u, x):
    a = 1
    e = 10
    iner = 1
    f = 10
    # Residual calculation
    r_d = a*e*sym.diff(sym.diff(u, x), x) + f
    return r_d


c1 = sym.Symbol('c1')
c2 = sym.Symbol('c2')
z = sym.Symbol('z')
c = []
# polynomial degree n-1
n = 10
u = 0
for i in range(n):
    c.append(sym.Symbol(f'c_{i}'))
    if i == 0:
        c[0] = 0
    u += (c[i]*z**i)
# u = u.subs({ : 0})
R_d = residual_calc(u, z)
print(R_d)
k = 0
soln_list = []
for i in range(1, n):
    k += i*c[i]
soln_list.append(k)
for i in range(n, 2, -1):
    temp_soln = R_d.subs({z: 1/i})
    soln_list.append(temp_soln)
soln = sym.solve(soln_list)
print(soln, c)
u_final = 0
for i in range(n):
    if i != 0:
        c[i] = soln[c[i]]
        u_final += c[i]*z**i
    else:
        u_final += c[i]*z**i
print(u_final)
# u = c1*z**3 + c2*z**2
# R_d = residual_calc(u, z)
# soln = sym.solve([R_d.subs({z: 0.5}), R_d.subs({z: 0.75})], dict= True)
# print(soln)
# c1 = soln[0][c1]
# c2 = soln[0][c2]
# u = c1*z**3 + c2*z**2
# print(u)
k = np.arange(0, 1, 0.1)
y = []
for num in k:
    u_temp = 0
    for j in range(n):
        if j != 0:
            u_temp += c[j]*num**j
        else:
            u_temp += c[j]*num**j
    y.append(u_temp)
plt.plot(k, y)
plt.show()
# print(residual_calc(u, z).subs({z: 10}).evalf())

