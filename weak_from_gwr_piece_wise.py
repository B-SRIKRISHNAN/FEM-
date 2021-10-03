
import sympy as sym
import math as mt
import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt


def residual_calc(u, x):
    # The area of cross_section
    a = 3E-5
    # Young's modulus
    e = 200E9
    # Area moment of inertia
    iner = 1
    # Residual calculation(done with integration of G.DE because we use weak form)
    r_d = a*e*sym.diff(u, x)
    return r_d


z = sym.Symbol('z')
L = 1
n = 4
ele_len = L / n
u = []
f_reac = []

# internal distributed forces( from G.D.E)

f_net = np.empty(shape=(n+1, 1), dtype=float)
for i in range(n+1):
    f_net[i] = 0

# Applied external forces

for i in range(n+1):
    if i == n/2+1:
        f_reac.append(100)
    else:
        f_reac.append(0)
print(f_reac)

for i in range(n+1):
    u.append(sym.Symbol(f'u{i+1}'))
print(u)
# u = 0
w1 = 1 - z / ele_len
w2 = z / ele_len
# Shape functions
N = [w1, w2]
node_list = []

for i in range(n):
    u_e = [u[i], u[i+1]]
    node_list.append(u_e)

R_d = []

for i in range(n):
    # The displacement interpolated within each element by shape functions
    ele_disp = np.dot(N, node_list[i])
    R_d.append(residual_calc(ele_disp, z))
print(R_d)

total_stiff = np.empty(shape=(n+1, n+1), dtype=float)
for i in range(n+1):
    for j in range(n+1):
        total_stiff[i][j] = 0

# distributed force
f = 100
# To get total_stiffness matrix of structure

for i in range(n):
    # Solving integral term in weak form equation for 1st weight function

    soln_temp = sym.integrate(sym.diff(N[0]) * R_d[i], (z, 0, ele_len))

    # Getting the coefficients of nodal displacements and applying them to total stiffness matrix

    c = sym.lambdify(node_list[i], soln_temp)
    total_stiff[i][i] += c(1, 0)
    total_stiff[i][i+1] += c(0, 1)

    # This gives us distributed force in nodes

    c_temp = sym.integrate(N[0] * (f if i < n/2 else 0), (z, 0, ele_len)).evalf()
    f_net[i] = f_net[i] + c_temp

    # Solving integral term in weak form equation for 2nd weight function

    soln_temp = sym.integrate(sym.diff(N[1]) * R_d[i], (z, 0, ele_len))

    # Getting the coefficients of nodal displacements and applying them to total stiffness matrix

    c = sym.lambdify(node_list[i], soln_temp)
    total_stiff[i+1][i] += c(1, 0)
    total_stiff[i+1][i+1] += c(0, 1)

    #  This gives us distributed force in nodes

    c_temp = sym.integrate(N[1] * (f if i < n/2 else 0), (z, 0, ele_len)).evalf()
    f_net[i+1] = f_net[i+1] + c_temp


print(f_net)
print(total_stiff)
k = np.empty(shape=(n+1, n+1))
for i in range(n+1):
    for j in range(n+1):
        k[i][j] = 0

# Applying natural boundary condition

for i in range(n+1):
    if i == 0:
        u[i] = 0

total_stiff_ad = total_stiff
total_stiff_pd = []
disp_ad = u
force_ad = []
for i in range(n+1):
    force_ad.append(float(f_net[i]) + float(f_reac[i]))
print(force_ad)

# getting the total stiffness and forces for active and prescribed dofs
for i in range(n+1):

    if u[i] == 0:
        total_stiff_ad = np.delete(total_stiff_ad, i, 0)
        total_stiff_ad = np.delete(total_stiff_ad, i, 1)
        disp_ad = np.delete(disp_ad, i, 0)
        force_ad = np.delete(force_ad, i, 0)
        total_stiff_pd.append(total_stiff[i])

print(total_stiff_ad)
print(force_ad)
total_stiff_ad_inv = np.linalg.inv(total_stiff_ad)
disp_ad_soln = np.dot(total_stiff_ad_inv, force_ad)
# print(disp_ad_soln)

for i in range(n+1):
    if i != 0:
        u[i] = float(disp_ad_soln[i-1])
# final nodal displacements
print(u)
force_pd = np.dot(total_stiff_pd, u) - float(f_net[0])

# reaction force at initial node
print(force_pd)

# Plotting the nodal displacements with position along bar
k = np.arange(0, L, 0.001)
disp = []
for num in k:
    for node in range(0, n):
        if num == node*ele_len:
            disp.append(u[node])
        elif node*ele_len < num < (node + 1)*ele_len:
            disp_temp = np.dot([1 - (num - node * ele_len) / ele_len, (num - node * ele_len) / ele_len], [u[node], u[node + 1]])
            disp.append(disp_temp)
        elif num == L:
            disp.append((u[n+1]))

plt.title('Trial_function - weak_form')
plt.xlabel('length')
plt.ylabel('Lateral displacement')
plt.plot(k, disp)
plt.show()

