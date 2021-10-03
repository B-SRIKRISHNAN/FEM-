import sympy as sym
import math as mt
import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt
G = 1
rho = 1
A = 1
E = 1
 # distributed force

x = sym.Symbol('x')
u1, u2, u3 = sym.symbols('u1, u2, u3')

num_ele = 1
num_nodes = 2*num_ele + 1
tot_len = 1
ele_len = tot_len/num_ele
Q = G*rho*A

# u1, u3, u2
shape_funcs = [1-3*x/ele_len + 2*(x**2)/ele_len**2, 4*x/ele_len - 4*(x**2)/ele_len**2, -x/ele_len + 2*(x**2)/ele_len**2]

f_dist = np.empty(shape=(num_nodes, 1), dtype=float)
for i in range(num_nodes):
    f_dist[i] = 0

f_nodal = np.empty(shape=(num_nodes, 1), dtype=float)
for i in range(num_nodes):
    f_nodal[i] = 0

disp_nodal = np.empty(shape=(num_nodes, 1), dtype=float)
for i in range(num_nodes):
    if i == 0:
        disp_nodal[i] = 0
    else:
        disp_nodal[i] = None


diff_shape_func = []


for func in shape_funcs:
    diff_shape_func.append(sym.diff(func, x))
# print(diff_shape_func)
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

btb = A*E*np.dot(bt, b)
# print(btb)
ele_stiff = []
for i in range(3):
    temp = []
    for j in range(3):
        temp.append(sym.integrate(btb[i][j], (x, 0, ele_len)))
    ele_stiff.append(temp)
# print(ele_stiff)
total_stiff = np.empty(shape=(num_nodes, num_nodes), dtype=float)
for i in range(num_nodes):
    for j in range(num_nodes):
        total_stiff[i][j] = 0
for i in range(num_ele):
    total_stiff[2*i][2*i] += ele_stiff[0][0]
    total_stiff[2*i][2*i+1] += ele_stiff[0][1]
    total_stiff[2*i][2*i+2] += ele_stiff[0][2]
    total_stiff[2*i+1][2*i] += ele_stiff[1][0]
    total_stiff[2*i+1][2*i+1] += ele_stiff[1][1]
    total_stiff[2*i+1][2*i+2] += ele_stiff[1][2]
    total_stiff[2*i+2][2*i] += ele_stiff[2][0]
    total_stiff[2*i+2][2*i+1] += ele_stiff[2][1]
    total_stiff[2*i+2][2*i+2] += ele_stiff[2][2]
    f_dist[2*i] = f_dist[2*i] + sym.integrate(shape_funcs[0]*Q, (x, 0, ele_len))
    f_dist[2*i+1] = f_dist[2*i+1] + sym.integrate(shape_funcs[1]*Q, (x, 0, ele_len))
    f_dist[2*i+2] = f_dist[2*i+2] + sym.integrate(shape_funcs[2]*Q, (x, 0, ele_len))
f_tot = f_nodal + f_dist
print(total_stiff)
# print(f_tot)

total_stiff_ad = total_stiff
total_stiff_pd = []
disp_ad = disp_nodal
force_ad = f_tot

for i in range(num_nodes):

    if disp_nodal[i] == 0:
        total_stiff_ad = np.delete(total_stiff_ad, i, 0)
        total_stiff_ad = np.delete(total_stiff_ad, i, 1)
        disp_ad = np.delete(disp_ad, i, 0)
        force_ad = np.delete(force_ad, i, 0)
        total_stiff_pd.append(total_stiff[i])

total_stiff_ad_inv = np.linalg.inv(total_stiff_ad)
disp_ad_soln = np.dot(total_stiff_ad_inv, force_ad)

for i in range(num_nodes):
    if i != 0:
        disp_nodal[i] = float(disp_ad_soln[i-1])
# final nodal displacements
print(disp_nodal)
force_pd = np.dot(total_stiff_pd, disp_nodal) - float(f_dist[0])

# reaction force at initial node
print(force_pd)
