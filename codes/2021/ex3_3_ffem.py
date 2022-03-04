# contoh 3.3 buku FFEM by Hutton
# 3D bar

import numpy as np
import calfem.core as cfc

Edof = np.array([
    [1, 2, 3, 10, 11, 12], # koord elemen 1
    [4, 5, 6, 10, 11, 12], # koord elemen 2
    [7, 8, 9, 10, 11, 12]  # koord elemen 3
])

coord = np.array([
    [0, 0, 30],
    [0, 0, -30], 
    [0, -30, 0], 
    [40, 0, 0]
])

Dof = np.array([
    [1, 2, 3], 
    [4, 5, 6],
    [7, 8, 9],
    [10, 11, 12]
])

ex, ey, ez = cfc.coordxtr(Edof, coord, Dof)

ep = np.array([3e5, 1])

# matriks kekakuan, lokal --> global
Ke1 = cfc.bar3e(ex[0], ey[0], ez[0], ep)
Ke1 = cfc.bar3e(ex[0], ey[0], ez[0], ep)