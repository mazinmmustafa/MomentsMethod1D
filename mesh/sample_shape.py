import numpy as np
import mesh
# 
lambda_ = 0.5
lc = 0.1
# 
lines_list = []
segments_list = []
# Lines
v1 = mesh.Vertex(0.0, 0.0, -0.5)
v2 = mesh.Vertex(0.0, 0.0, +0.5)
v3 = mesh.Vertex(+0.5, 0.0, +0.5)
v4 = mesh.Vertex(-0.5, 0.0, +0.5)
lines_list.append(mesh.Line(v1, v2, True, 1))
lines_list.append(mesh.Line(v3, v2, False, 0))
lines_list.append(mesh.Line(v4, v2, False, 0))
# Circle
radius = 0.5
N = int(np.ceil(2.0*np.pi*(radius/lambda_)/lc))
if N<6:
    N = 6
    pass 
d_phi = 2.0*np.pi/N
for i in range(N):
    phi = (i+0)*d_phi
    v1 = mesh.Vertex(radius*np.cos(phi), radius*np.sin(phi), +0.5)
    phi = (i+1)*d_phi
    v2 = mesh.Vertex(radius*np.cos(phi), radius*np.sin(phi), +0.5)
    if i==0 or i==N-1:
        segment = mesh.Segment(v1, v2, True, 2)
        pass
    else:
        segment = mesh.Segment(v1, v2, False, 0)
        pass
    segments_list.append(segment)
    continue
# 
shape = mesh.Shape(lambda_, lc)
for line in lines_list:
    shape.add_line(line)
    continue
for segment in segments_list:
    shape.add_segment(segment)
    continue
shape.mesh()

