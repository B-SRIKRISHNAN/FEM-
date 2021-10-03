# Modified from truss_gen, This generalises x and y inputs. You can give separate
# x and y displacements for a node or just give one and leave other as an active dof

import numpy as np
import math as mt
num_ele = int(input('Enter number of elements'))
node_ele_rel = []   # Element nodal relation
stiff = []          # To collect the element stiffness
disp = []           # storing the nodal displacements
force = []          # Storing the nodal forces
angle = []
# force_nh = []       # Non homogenous force matrix
# To get the stiffness and nodes of each element
for i in range(num_ele):

    youngs_mod = 100E9
    length = float(input(f'Enter length of element {i+1}\t'))
    area_cross = 1E-4
    stiffness = youngs_mod*area_cross/length
    stiff.append(stiffness)
    angle.append(float(input(f'Enter angle of element {i+1} in degrees\t')))
    print(f'Enter the nodes of element {i+1}\t')
    node1 = int(input('first node\t'))
    node2 = int(input('second node\t'))
    node_ele_rel.append([node1, node2])

trans_mat = []
for k in range(num_ele):
    temp = np.empty(shape=(2, 4), dtype=float)
    temp[0][0] = round(mt.cos(mt.radians(angle[k])),1)
    temp[0][1] = round(mt.sin(mt.radians(angle[k])),1)
    temp[0][2] = 0
    temp[0][3] = 0
    temp[1][0] = 0
    temp[1][1] = 0
    temp[1][2] = round(mt.cos(mt.radians(angle[k])), 1)
    temp[1][3] = round(mt.sin(mt.radians(angle[k])), 1)

    trans_mat.append(temp)

num_nodes = 0

# To get the number of nodes

for i in range(num_ele):
    for j in range(2):
        if node_ele_rel[i][j] > num_nodes:
            num_nodes = node_ele_rel[i][j]
# print(num_nodes)

total_stiff = np.empty(shape=(num_nodes*2, num_nodes*2), dtype=float)
for i in range(num_nodes*2):
    for j in range(num_nodes*2):
        total_stiff[i][j] = 0

# To get the total stiffness matrix
ele_stiff = []
for k in range(num_ele):
    temp = np.empty(shape=(2,2), dtype=float)
    temp[0][0] = stiff[k]
    temp[0][1] = -stiff[k]
    temp[1][0] = -stiff[k]
    temp[1][1] = stiff[k]
    trans_mat_trans = trans_mat[k].transpose()
    print(trans_mat_trans)
    temp2 = np.dot(trans_mat_trans, temp)

    final_stiff_mat = np.dot(temp2, trans_mat[k])

    ele_stiff.append(final_stiff_mat)
ele_stiff = np.array(ele_stiff)
# print('The element stiffness matrices are\n', ele_stiff)
# print(ele_stiff[0][0][1])
for k in range(num_ele):
    for i in range(num_nodes):
        for j in range(num_nodes):
            if [i+1, j+1] == node_ele_rel[k]:
                # stiffness = stiff[node_ele_rel.index([i+1, j+1])]
                total_stiff[2*i][2*i] += ele_stiff[k][0][0]
                total_stiff[2*i][2*i+1] += ele_stiff[k][0][1]
                total_stiff[2*i+1][2*i] += ele_stiff[k][1][0]
                total_stiff[2*i+1][2*i+1] += ele_stiff[k][1][1]
                total_stiff[2*i][2*j] += ele_stiff[k][0][2]
                total_stiff[2*i][2*j+1] += ele_stiff[k][0][3]
                total_stiff[2*i+1][2*j] += ele_stiff[k][1][2]
                total_stiff[2*i+1][2*j+1] += ele_stiff[k][1][3]
                total_stiff[2*j][2*i] += ele_stiff[k][2][0]
                total_stiff[2*j][2*i+1] += ele_stiff[k][2][1]
                total_stiff[2*j+1][2*i] += ele_stiff[k][3][0]
                total_stiff[2*j+1][2*i+1] += ele_stiff[k][3][1]
                total_stiff[2*j][2*j] += ele_stiff[k][2][2]
                total_stiff[2*j][2*j+1] += ele_stiff[k][2][3]
                total_stiff[2*j+1][2*j] += ele_stiff[k][3][2]
                total_stiff[2*j+1][2*j+1] += ele_stiff[k][3][3]
print('The total stiffness matrix is \n', total_stiff)

# Getting the nodal displacements

for i in range(num_nodes):
    choice = input(f'Enter y if x displacement is given for node {i+1}\t')
    if choice == 'y':
        disp_x = float(input(f'enter x displacement of node {i+1}\t'))
    else:
        disp_x = None
    choice = input(f'Enter y if y displacement is given for node{i+1}\t')
    if choice == 'y':
        disp_y = float(input(f'Enter y displacement of node {i+1}\t'))
    else:
        disp_y = None
    disp.append([disp_x, disp_y])


disp_fin = []
for i in range(num_nodes):
    disp_fin.append(disp[i][0])
    disp_fin.append(disp[i][1])
disp = disp_fin


total_stiff_ad = total_stiff
total_stiff_ad_pd = total_stiff
total_stiff_pd = []
# Getting the nodal forces

for i in range(num_nodes):

    choice = input(f'Enter y if x - force exists for node{i+1}\t')
    if choice == 'y':
        force_x = float(input(f'enter x force of node {i+1}\t'))
    else:
        force_x = None
    choice = input(f'Enter y if y - force exists for node{i+1}\t')
    if choice == 'y':
        force_y = float(input(f'enter y force of node {i+1}\t'))
    else:
        force_y = None

    force.append([force_x, force_y])

force_fin = []
for i in range(num_nodes):
    force_fin.append(force[i][0])
    force_fin.append(force[i][1])
force = force_fin
# print('Input displacements\n', disp)
# print('Input force\n', force)

num_pd = 0
disp_pd = []
num_ad = 0
force_ad = np.array(force)
for i in range(num_nodes*2):
    if disp[i] is None:
        total_stiff_ad_pd = np.delete(total_stiff_ad_pd, (i - num_ad), 1)
        num_ad += 1
    else:
        disp_pd.append(disp[i])

        total_stiff_pd.append(total_stiff[i])

        total_stiff_ad = np.delete(total_stiff_ad, (i - num_pd), 0)

        total_stiff_ad = np.delete(total_stiff_ad, (i - num_pd), 1)

        total_stiff_ad_pd = np.delete(total_stiff_ad_pd, (i - num_pd), 0)

        force_ad = np.delete(force_ad, (i - num_pd), 0)

        num_pd += 1

# print('Total_stiff active dof', total_stiff_ad)
# force_ad = np.array(force_ad)
force_ad_nh = force_ad - np.dot(total_stiff_ad_pd, disp_pd)
# print(force_ad_nh)
total_stiff_ad_inv = np.linalg.inv(total_stiff_ad)
disp_ad = np.dot(total_stiff_ad_inv, force_ad_nh)

num_pd = 0
for i in range(num_nodes*2):
    if disp[i] is None:
        disp[i] = disp_ad[(i - num_pd)]

    else:
        num_pd += 1


disp = np.array(disp)
force_pd = np.dot(total_stiff_pd, disp)
n = 0
for i in range(num_nodes*2):
    if force[i] is None:
        force[i] = force_pd[n]
        n += 1


print('Nodal forces are\n', force)
print('Nodal displacements are\n', disp)
