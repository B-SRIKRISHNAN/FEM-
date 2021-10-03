# This is a general code for spring stiffness in 1-d

import numpy as np
num_ele = int(input('Enter number of elements'))
node_ele_rel = []   # Element nodal relation
stiff = []          # To collect the element stiffness
disp = []           # storing the nodal displacements
force = []          # Storing the nodal forces

# force_nh = []       # Non homogenous force matrix
# To get the stiffness and nodes of each element

for i in range(num_ele):
    stiff.append(float(input(f'Enter elemental stiffness for element {i+1}')))
    print(f'Enter the nodes of element {i+1}')
    node1 = int(input('first node'))
    node2 = int(input('second node'))
    node_ele_rel.append([node1, node2])

num_nodes = 0

# To get the number of nodes

for i in range(num_ele):
    for j in range(2):
        if node_ele_rel[i][j] > num_nodes:
            num_nodes = node_ele_rel[i][j]
# print(num_nodes)

total_stiff = np.empty(shape=(num_nodes, num_nodes), dtype=float)
for i in range(num_nodes):
    for j in range(num_nodes):
        total_stiff[i][j] = 0

# To get the total stiffness matrix
ele_stiff = []
for k in range(num_ele):
    temp = np.empty(shape=(2,2), dtype=float)
    temp[0][0] = stiff[k]
    temp[0][1] = -stiff[k]
    temp[1][0] = -stiff[k]
    temp[1][1] = stiff[k]
    ele_stiff.append(temp)
ele_stiff = np.array(ele_stiff)

print('The element stiffness matrices are\n', ele_stiff)

# print(ele_stiff[0][0][1])

for k in range(num_ele):
    for i in range(num_nodes):
        for j in range(num_nodes):
            if [i+1, j+1] == node_ele_rel[k]:
                # stiffness = stiff[node_ele_rel.index([i+1, j+1])]
                total_stiff[i][i] += ele_stiff[k][0][0]
                total_stiff[i][j] += ele_stiff[k][0][1]
                total_stiff[j][i] += ele_stiff[k][1][0]
                total_stiff[j][j] += ele_stiff[k][1][1]

print('The total stiffness matrix is \n', total_stiff)

# Getting the nodal displacements

for i in range(num_nodes):
    choice = input(f'Enter y if there is nodal displacement for node {i+1}, else just press enter')
    if choice == 'y':
        disp.append(float(input(f'Enter nodal displacement for node {i+1}')))
    else:
        disp.append(None)
total_stiff_ad = total_stiff
total_stiff_ad_pd = total_stiff
total_stiff_pd = []
# Getting the nodal forces

for i in range(num_nodes):

    if disp[i] is None:
        force.append(float(input(f'Enter nodal force for node {i+1}')))
    else:
        force.append(None)


num_pd = 0
disp_pd = []
num_ad = 0
force_ad = force
for i in range(num_nodes):
    if disp[i] is None:
        total_stiff_ad_pd = np.delete(total_stiff_ad_pd, i - num_ad, 1)
        num_ad += 1
    else:
        disp_pd.append(disp[i])
        total_stiff_pd.append(total_stiff[i])
        total_stiff_ad = np.delete(total_stiff_ad, i - num_pd, 0)
        total_stiff_ad = np.delete(total_stiff_ad, i - num_pd, 1)
        total_stiff_ad_pd = np.delete(total_stiff_ad_pd, i - num_pd, 0)
        force_ad = np.delete(force_ad, i - num_pd, 0)
        num_pd += 1


force_ad = np.array(force_ad)
force_ad_nh = force_ad - np.dot(total_stiff_ad_pd, disp_pd)
total_stiff_ad_inv = np.linalg.inv(total_stiff_ad)
disp_ad = np.dot(total_stiff_ad_inv, force_ad_nh)

num_pd = 0
for i in range(num_nodes):
    if disp[i] is None:
        disp[i] = disp_ad[i - num_pd]

    else:
        num_pd += 1


disp = np.array(disp)
force_pd = np.dot(total_stiff_pd, disp)
n = 0
for i in range(num_nodes):
    if force[i] is None:
        force[i] = force_pd[n]
        n += 1


print('Nodal forces are\n', force)
print('Nodal displacements are\n', disp)
