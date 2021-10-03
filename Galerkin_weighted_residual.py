
import sympy as sym
import math as mt
import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt


def residual_calc(u, x):
    a = 1
    e = 10
    iner = 1
    f = 10
    # Residual calculation
    r_d = iner*e*sym.diff(sym.diff(sym.diff(sym.diff(u, x), x), x),  x) - f
    return sym.sin(np.pi*x)*r_d


z = sym.Symbol('z')
c1 = sym.Symbol('c1')
# c2 = sym.Symbol('c2')
# c0 = sym.Symbol('c0')
# polynomial degree n-1
# n = 1
# u = 0

# Trial solution
u = c1*sym.sin(np.pi*z)


R_d = residual_calc(u, z)
print(R_d)

soln_list = sym.integrate(R_d, (z, 0, 1))
print(soln_list.evalf())
soln = sym.solve(soln_list.evalf())
print(soln[0])
c1 = soln[0]
# u_final = c1*sym.sin(np.pi*z)

#
k = np.arange(0, 1, 0.01)
y = []
for num in k:
    u_temp = c1*sym.sin(np.pi*num)
    y.append(u_temp)
plt.title('Galerkin Weighted Residual bar')
plt.xlabel('length')
plt.ylabel('Lateral displacement')
plt.plot(k, y)
plt.show()

