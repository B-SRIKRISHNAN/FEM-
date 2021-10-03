import sympy as sym
import math as mt
import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt
G = 1
rho = 1
A = 1
E = 1
I = 1  # second moment of area
# distributed force
# z = sym.Symbol('z')
x = sym.Symbol('x')
u1, u2, u3 = sym.symbols('u1, u2, u3')

num_ele = 5
num_nodes = num_ele + 1
node_dof = 2
tot_len = 5
ele_len = tot_len/num_ele
Q = -12E3

# v1, theta1, v2, theta2
shape_funcs = [1-(3*x**2)/(ele_len**2) + 2*(x**3)/ele_len**3, x - (2*x**2/ele_len) + x**3/ele_len**3, 3*x**2/ele_len**2
               - 2*(x**3)/ele_len**3, -x**2/ele_len + (x**3)/ele_len**2]

f_dist = np.empty(shape=(num_nodes*node_dof, 1), dtype=float)
for i in range(num_nodes*node_dof):
    f_dist[i] = 0

f_nodal = np.empty(shape=(num_nodes*node_dof, 1), dtype=float)
for i in range(num_nodes*node_dof):
    if i == 2:
        f_nodal[i] = -1000
    else:
        f_nodal[i] = 0

disp_nodal = np.empty(shape=(num_nodes*node_dof, 1), dtype=float)
for i in range(num_nodes*node_dof):
    if i == 0 or i == 1 or i == 4 or i == 6 or i == 10:
        disp_nodal[i] = 0
    else:
        disp_nodal[i] = None


diff_shape_func = []


for func in shape_funcs:
    diff_shape_func.append(sym.diff(sym.diff(func, x), x))
print(diff_shape_func)
# bt = []
bt = np.empty(shape=(len(shape_funcs), 1), dtype=object)
for i in range(len(shape_funcs)):
    bt[i][0] = (diff_shape_func[i])
# print(bt)
# bt = np.array(bt)
# bt = np.transpose(bt)
b = np.empty(shape=(1, len(shape_funcs)), dtype=object)
for i in range(len(shape_funcs)):
    b[0][i] = (diff_shape_func[i])
# print(b)
#
btb = I*E*np.dot(bt, b)
# print(btb)
ele_stiff = []
for i in range(2*node_dof):
    temp = []
    for j in range(2*node_dof):
        temp.append(sym.integrate(btb[i][j], (x, 0, ele_len)))
    ele_stiff.append(temp)
# print(ele_stiff)
total_stiff = np.empty(shape=(num_nodes*node_dof, num_nodes*node_dof), dtype=float)
for i in range(num_nodes*node_dof):
    for j in range(num_nodes*node_dof):
        total_stiff[i][j] = 0
for i in range(num_ele):
    for k in range(4):
        for l in range(4):
            total_stiff[2*i+k][2*i+l] += ele_stiff[0+k][0+l]
        f_dist[2*i+k] = f_dist[2*i+k] + sym.integrate(shape_funcs[0+k]*(Q if (i >= 3) else 0), (x, 0, ele_len))


# f_dist[2*i+1] = f_dist[2*i+1] + sym.integrate(shape_funcs[1]*Q, (x, 0, ele_len))
    # f_dist[2*i+2] = f_dist[2*i+2] + sym.integrate(shape_funcs[2]*Q, (x, 0, ele_len))
f_tot = f_nodal + f_dist
# print(total_stiff)
print(f_tot)
#
total_stiff_ad = total_stiff
total_stiff_pd = []
disp_ad = disp_nodal
force_ad = f_tot
#
n = 0
for i in range(num_nodes*node_dof):

    if disp_nodal[i] == 0:
        total_stiff_ad = np.delete(total_stiff_ad, i-n, 0)
        total_stiff_ad = np.delete(total_stiff_ad, i-n, 1)
        disp_ad = np.delete(disp_ad, i-n, 0)
        force_ad = np.delete(force_ad, i-n, 0)
        total_stiff_pd.append(total_stiff[i])
        n += 1
# print('total_stiff_ad', total_stiff_ad)
total_stiff_ad_inv = np.linalg.inv(total_stiff_ad)
disp_ad_soln = np.dot(total_stiff_ad_inv, force_ad)
n2 = 0
for i in range(num_nodes*node_dof):
    if disp_nodal[i] != 0:
        disp_nodal[i] = float(disp_ad_soln[i-n2])
    else:
        n2 += 1
# final nodal displacements
print(disp_nodal)
print(f_dist)
force_pd = np.array(np.dot(total_stiff_pd, disp_nodal)) \
        - np.array([f_dist[0], f_dist[1], f_dist[4], f_dist[6], f_dist[10]])
# reaction force at initial node
print(np.array(force_pd))
