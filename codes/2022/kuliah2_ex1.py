import numpy as np
import calfem.core as cfc

Edof = np.array([
    [1, 2],      # element 1 between node 1 and 2
    [2, 3],      # element 2 between node 2 and 3
    [2, 3]       # element 3 between node 2 and 3
])

K = np.zeros((3,3))
f = np.zeros((3,1))

# Element stiffness matrices  
k = 1500.0
ep1 = k                  
ep2 = 2.*k
Ke1 = cfc.spring1e(ep1)
Ke2 = cfc.spring1e(ep2)

# Masukkan nilai Ke untuk masing-masing elemen
cfc.assem(Edof[0,:], K, Ke2)        # element 1
cfc.assem(Edof[1,:], K, Ke1)        # element 2
cfc.assem(Edof[2,:], K, Ke2)        # element 3

# definisikan gaya
f[1] = 100

# F = k*x = K * U = K * d

# set up boundary conditions
bc = np.array([1,3])  # node 1 and 3 are constrained

# solve for displacement
# a is nodal displacement, r is reaction force
a, r = cfc.solveq(K, f, bc)

# for finding element displacements
ed1 = cfc.extractEldisp(Edof[0,:], a)
ed2 = cfc.extractEldisp(Edof[1,:], a)
ed3 = cfc.extractEldisp(Edof[2,:], a)

# Untuk mencari gaya pada tiap elemen/pegas
es1 = cfc.spring1s(ep2, ed1)
es2 = cfc.spring1s(ep1, ed2)
es3 = cfc.spring1s(ep2, ed3)


# gaya normal untuk tiap elemen
print(f"N1 = {es1}")
print(f"N2 = {es2}")
print(f"N3 = {es3}")