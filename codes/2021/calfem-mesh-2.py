# calfem tutorial demo: complete mesh with solver
# from calfem documentation: 

import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv
import calfem.utils as cfu

import numpy as np

# Describing simple rectangular beam
l = 5.0             # length
h = 1.0             # height 
t = 0.2             # thickness

# define constants for marker
left_support = 10
right_support = 20
top_line = 30

# Creating geometry
g = cfg.Geometry()

g.point([0.0, 0.0], marker = left_support) # point 0
g.point([l, 0.0], marker = right_support) # point 1
g.point([l, h]) # point 2
g.point([0.0, h]) # point 3

g.spline([0, 1]) # line 0
g.spline([1, 2]) # line 1
g.spline([2, 3], marker = top_line) # line 2
g.spline([3, 0]) # line 3

g.surface([0, 1, 2, 3])

# Creating a mesh
mesh = cfm.GmshMesh(g)

mesh.elType = 2          # Degrees of freedom per node.
mesh.dofsPerNode = 1     # Factor that changes element sizes.
mesh.elSizeFactor = 0.15 # Element size Factor

coords, edof, dofs, bdofs, elementmarkers = mesh.create()

# Draw the mesh.
cfv.drawMesh(
    coords=coords,
    edof=edof,
    dofs_per_node=mesh.dofsPerNode,
    el_type=mesh.elType,
    filled=True,
    title="Example 01"
    )