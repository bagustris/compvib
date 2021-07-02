#  contoh 3.7 buku FFEM by Hutton

import numpy as np
import calfem.core as cfc

# definisikan Edof
Edof = np.array([
    [1, 2, 5, 6],     # elemen 1
    [1, 2, 7, 8],     # elemen 2
    [3, 4, 7, 8],     # elemen 3
    [5, 6, 7, 8],     # elemen 4
    [5, 6, 9, 10],    # elemen 5
    [7, 8, 9, 10],    # elemen 6  
    [7, 8, 11, 12],   # elemen 7
    [9, 10, 11, 12]   # elemen 8 
])

Dof= np.array([
    [1, 2], 
    [3, 4],
    [5, 6], 
    [7, 8], 
    [9, 10], 
    [11, 12] 
])

coord = np.array([
    [0, 0], 
    [0, 40], 
    [40, 0], 
    [40, 40], 
    [80, 0], 
    [80, 40]
])


ex, ey = cfc.coordxtr(Edof, coord, Dof)


K = np.zeros((12, 12))
f = np.zeros((12, 1))
f[5] = -2000
f[8] = 2000
f[10] = 4000
f[11] = 6000

ep = np.array([10e6, 1.5])

# hitung K
Ke1 = cfc.bar2e(ex[0], ey[0], ep)
Ke2 = cfc.bar2e(ex[1], ey[1], ep)
Ke3 = cfc.bar2e(ex[2], ey[2], ep)
Ke4 = cfc.bar2e(ex[3], ey[3], ep)
Ke5 = cfc.bar2e(ex[4], ey[4], ep)
Ke6 = cfc.bar2e(ex[5], ey[5], ep)
Ke7 = cfc.bar2e(ex[6], ey[6], ep)
Ke8 = cfc.bar2e(ex[7], ey[7], ep)

cfc.assem(Edof[0], K, Ke1)
cfc.assem(Edof[1], K, Ke2)
cfc.assem(Edof[2], K, Ke3)
cfc.assem(Edof[3], K, Ke4)
cfc.assem(Edof[4], K, Ke5)
cfc.assem(Edof[5], K, Ke6)
cfc.assem(Edof[6], K, Ke7)
cfc.assem(Edof[7], K, Ke8)


bc = np.array([1, 2, 3, 4])
a, r = cfc.solveq(K, f, bc)

print("U= ", a)
print ("Gaya reaksi=", r)

ed = cfc.extractEldisp(Edof,a);
N = np.zeros([Edof.shape[0]])
s = np.zeros([Edof.shape[0]])
st = np.zeros([Edof.shape[0]])

print("Element forces:")

# hitung
i = 0
for elx, ely, eld in zip(ex, ey, ed):    
    N[i] = cfc.bar2s(elx,ely,ep,eld)
    s[i] = N[i]/1.5
    st[i] = s[i]/ep[0]
    print("stress_%d = %f" % (i+1,s[i]))
    print("strain_%d = %f" % (i+1,st[i]))
    # print("N%d = %g" % (i+1,N[i]))
    i+=1