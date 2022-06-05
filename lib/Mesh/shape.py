import Mesh1D as msh
# Definitions
basisFileName = 'Mesh/Basis.dat'
msehFileName = 'Mesh/Mesh.dat'
delta =  4.83870967741935E-03
# Define vertices
v1 = msh.Vertex(  0.00000000000000E+00,  0.00000000000000E+00, -3.00000000000000E-02)
v2 = msh.Vertex(  0.00000000000000E+00,  0.00000000000000E+00,  3.00000000000000E-02)
v3 = msh.Vertex(  1.00000000000000E-01,  0.00000000000000E+00,  3.00000000000000E-02)
v4 = msh.Vertex(  1.00000000000000E-01,  1.00000000000000E-01,  3.00000000000000E-02)
v5 = msh.Vertex(  1.00000000000000E-01,  1.00000000000000E-01, -3.00000000000000E-02)
v6 = msh.Vertex(  1.00000000000000E-01,  0.00000000000000E+00, -3.00000000000000E-02)
# Define lines
Line1 = msh.Line(v1, v2, 1, delta)
Line2 = msh.Line(v2, v3, 0, delta)
Line3 = msh.Line(v3, v4, 0, delta)
Line4 = msh.Line(v4, v5, 0, delta)
Line5 = msh.Line(v5, v6, 0, delta)
Line6 = msh.Line(v6, v1, 0, delta)
Line7 = msh.Line(v3, v6, 0, delta)
# Define shape
Shape1 = msh.Shape([Line1, Line2, Line3, Line4, Line5, Line6, Line7])
# Create mesh
myMesh = msh.meshShape(Shape1, msehFileName, True)
# Create Basis
N_Basis = msh.createBasis(myMesh, basisFileName)
# Save No. of basis
with open('Mesh/BasisInfo.txt', 'w') as file:
	file.write('{}'.format(N_Basis))
	pass
