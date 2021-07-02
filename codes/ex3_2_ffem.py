# contoh 3.2 buku FFEM by Hutton

import numpy as np
import calfem.core as cfc

Edof = np.array([
    [1, 2, 5, 6],  # element 1
    [3, 4, 5, 6]   # element 2
])

# koordinat utk ex, ey, bisa juga manual spt di bawahnya
coord = np.array([
    [0, 0],
    [0, 40],  
    [40, 40]
])

Dof = np.array([
    [1, 2], 
    [3, 4],  
    [5, 6]
])
ex, ey = cfc.coordxtr(Edof, coord, Dof)

# bisa juga spt ini, per elemen
ex1 = np.array([0, 40])
ey1 = np.array([0, 40])
ex2 = np.array([0, 40])
ey2 = np.array([40, 40])

K = np.zeros((6, 6))
f = np.zeros((6, 1))
f[4] = 500
f[5] = 300

# modulus young dan luas penampang
ep = [1e7, 1.5]

# hitung K total
Ke1 = cfc.bar2e(ex1, ey1, ep)
Ke2 = cfc.bar2e(ex2, ey2, ep)

cfc.assem(Edof[0], K, Ke1)
cfc.assem(Edof[1], K, Ke2)

# kondisi pembatas
bc = np.array([1, 2, 3, 4])
a, r = cfc.solveq(K, f, bc)

# hitung element force pada tiap elemen
ed = cfc.extractEldisp(Edof,a);
N = np.zeros([Edof.shape[0]])

print("Element forces:")

i = 0
for elx, ely, eld in zip(ex, ey, ed):
    N[i] = cfc.bar2s(elx,ely,ep,eld);
    print("N%d = %g" % (i+1,N[i]))
    i+=1

# calculate stress = N[i]/A
stress1 = N[0]/1.5
stress2 = N[1]/1.5