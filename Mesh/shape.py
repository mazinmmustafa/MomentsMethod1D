import Mesh1D as msh
import math
# Definitions
basisFileName = 'Mesh/Basis.dat'
msehFileName = 'Mesh/Mesh.dat'
delta =  3.22580645161290E-02
# Define vertices
v1 = msh.Vertex(  0.00000000000000E+00,  0.00000000000000E+00,  4.35000000000000E-02)
v2 = msh.Vertex(  0.00000000000000E+00,  0.00000000000000E+00,  4.96125000000000E-01)
v3 = msh.Vertex(  0.00000000000000E+00,  0.00000000000000E+00,  6.47000000000000E-01)
v4 = msh.Vertex(  0.00000000000000E+00,  0.00000000000000E+00, -4.35000000000000E-02)
v5 = msh.Vertex(  0.00000000000000E+00,  0.00000000000000E+00, -4.96125000000000E-01)
v6 = msh.Vertex(  0.00000000000000E+00,  0.00000000000000E+00, -6.47000000000000E-01)
Lines = []
Lines.append(msh.Line(v4, v1, 1, delta/3.0))
Lines.append(msh.Line(v1, v2, 0, delta))
Lines.append(msh.Line(v2, v3, 0, delta))
Lines.append(msh.Line(v5, v4, 0, delta))
Lines.append(msh.Line(v6, v5, 0, delta))
v7 = msh.Vertex(  2.61323165591954E-01,  0.00000000000000E+00,  4.96125000000000E-01)
Lines.append(msh.Line(v1, v7, 0, delta))
Lines.append(msh.Line(v7, v3, 0, delta))
v8 = msh.Vertex(  2.61323165591954E-01,  0.00000000000000E+00, -4.96125000000000E-01)
Lines.append(msh.Line(v4, v8, 0, delta))
Lines.append(msh.Line(v8, v6, 0, delta))
v9 = msh.Vertex(  1.30661582795977E-01,  2.26312500000000E-01,  4.96125000000000E-01)
Lines.append(msh.Line(v1, v9, 0, delta))
Lines.append(msh.Line(v9, v3, 0, delta))
v10 = msh.Vertex(  1.30661582795977E-01,  2.26312500000000E-01, -4.96125000000000E-01)
Lines.append(msh.Line(v4, v10, 0, delta))
Lines.append(msh.Line(v10, v6, 0, delta))
v11 = msh.Vertex( -1.30661582795977E-01,  2.26312500000000E-01,  4.96125000000000E-01)
Lines.append(msh.Line(v1, v11, 0, delta))
Lines.append(msh.Line(v11, v3, 0, delta))
v12 = msh.Vertex( -1.30661582795977E-01,  2.26312500000000E-01, -4.96125000000000E-01)
Lines.append(msh.Line(v4, v12, 0, delta))
Lines.append(msh.Line(v12, v6, 0, delta))
v13 = msh.Vertex( -2.61323165591954E-01,  3.20018008984885E-17,  4.96125000000000E-01)
Lines.append(msh.Line(v1, v13, 0, delta))
Lines.append(msh.Line(v13, v3, 0, delta))
v14 = msh.Vertex( -2.61323165591954E-01,  3.20018008984885E-17, -4.96125000000000E-01)
Lines.append(msh.Line(v4, v14, 0, delta))
Lines.append(msh.Line(v14, v6, 0, delta))
v15 = msh.Vertex( -1.30661582795977E-01, -2.26312500000000E-01,  4.96125000000000E-01)
Lines.append(msh.Line(v1, v15, 0, delta))
Lines.append(msh.Line(v15, v3, 0, delta))
v16 = msh.Vertex( -1.30661582795977E-01, -2.26312500000000E-01, -4.96125000000000E-01)
Lines.append(msh.Line(v4, v16, 0, delta))
Lines.append(msh.Line(v16, v6, 0, delta))
v17 = msh.Vertex(  1.30661582795977E-01, -2.26312500000000E-01,  4.96125000000000E-01)
Lines.append(msh.Line(v1, v17, 0, delta))
Lines.append(msh.Line(v17, v3, 0, delta))
v18 = msh.Vertex(  1.30661582795977E-01, -2.26312500000000E-01, -4.96125000000000E-01)
Lines.append(msh.Line(v4, v18, 0, delta))
Lines.append(msh.Line(v18, v6, 0, delta))
# Define lines
# Define shape
Shape1 = msh.Shape(Lines)
# Create mesh
myMesh = msh.meshShape(Shape1, msehFileName, True)
# Create Basis
N_Basis = msh.createBasis(myMesh, basisFileName)
# Save No. of basis
with open('Mesh/BasisInfo.txt', 'w') as file:
	file.write('{}'.format(N_Basis))
	pass
