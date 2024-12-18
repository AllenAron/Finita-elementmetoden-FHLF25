import numpy as np
import calfem.core as cfc

L = 6
num_elements = 3
dx = L/num_elements
num_nodes = num_elements + 1

A = 10
k = 5
Q = 100

Edof = np.array([[i + 1, i + 2] for i in range(num_elements)])

K = np.zeros((num_nodes, num_nodes))
F = np.zeros((num_nodes, 1))
fe = (Q * dx/2) * np.array([1, 1])

for i in range(num_elements):
    ke = cfc.spring1e(k*A*dx)
    cfc.assem(Edof[i], K, ke)
    F[Edof[i] - 1, 0] += fe

bc = np.array([1])
bcVal = np.array([0])
F[num_nodes - 1] += -A*15

a, r = cfc.solveq(K, F, bc, bcVal)
print(a)