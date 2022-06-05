import math

# Class Vertex
class Vertex:
    def __init__(self, x_in=0, y_in=0, z_in=0):
        self.x = x_in
        self.y = y_in
        self.z = z_in
        pass
    def __repr__(self):
        return "({:17.14f}, ".format(self.x)\
               +"{:17.14f}, ".format(self.y)\
               +"{:17.14f})".format(self.z)
    def __add__(self, other):
        return Vertex(self.x+other.x, self.y+other.y, self.z+other.z)
    def __sub__(self, other):
        return Vertex(self.x-other.x, self.y-other.y, self.z-other.z)
    def scale(self, a):
        return Vertex(a*self.x, a*self.y, a*self.z)

# Class Line made of verticies
class Line:
    def __init__(self, v1_in=Vertex(), v2_in=Vertex(), port_in=0, delta_in=0):
        self.v1 = v1_in
        self.v2 = v2_in
        self.port = port_in
        self.delta = delta_in
        pass
    def __repr__(self):
        return "v1 = ({:17.14f}, ".format(self.v1.x)\
               +"{:17.14f}, ".format(self.v1.y)\
               +"{:17.14f})\n".format(self.v1.z)\
               +"v2 = ({:17.14f}, ".format(self.v2.x)\
               +"{:17.14f}, ".format(self.v2.y)\
               +"{:17.14f})".format(self.v2.z)
    def len(self):
       return math.sqrt(math.pow(self.v1.x-self.v2.x, 2)\
                       +math.pow(self.v1.y-self.v2.y, 2)\
                       +math.pow(self.v1.z-self.v2.z, 2))

# Class Shape made of lines
class Shape:
    def __init__(self, Lines_in=[]):
        self.Lines = []
        for line in Lines_in:
            self.Lines.append(line)
            continue
        pass
    def __repr__(self):
        str = ""
        count = 0
        for line in self.Lines:
            count+=1
            str = str+"Line {:d}\n".format(count)\
                   +"v1 = ({:17.14f}, ".format(line.v1.x)\
                   +"{:17.14f}, ".format(line.v1.y)\
                   +"{:17.14f})\n".format(line.v1.z)\
                   +"v2 = ({:17.14f}, ".format(line.v2.x)\
                   +"{:17.14f}, ".format(line.v2.y)\
                   +"{:17.14f})\n".format(line.v2.z)
            continue
        return str

# Create a mesh from shape
def meshShape(myShape, fileName):
    meshList = []
    for line in myShape.Lines:
        L = line.len()
        delta = line.delta
        assert(delta>0)
        if delta>L:
            delta = L
            pass
        n = L/delta
        n = math.ceil(n)
        if n%2 != 0:
            n+=1
            pass
        delta = L/n
        for i in range(0, n):
            alpha_i = (i+0)*delta/L
            node1 = Vertex()
            node1 = line.v1+(line.v2-line.v1).scale(alpha_i)
            alpha_i = (i+1)*delta/L
            node2 = Vertex()
            node2 = line.v1+(line.v2-line.v1).scale(alpha_i)
            if line.port!=0 and (i==(n/2-1) or i==(n/2)):
                newEntrie = [ node1.x, node1.y, node1.z, \
                              node2.x, node2.y, node2.z, line.port ]
                pass
            else:
                newEntrie = [ node1.x, node1.y, node1.z, \
                              node2.x, node2.y, node2.z, 0 ]
                pass
            meshList.append(newEntrie)
            continue
        continue
    n = len(meshList)
    # print("{} elements were found!".format(n))
    with open(fileName, "w") as file:
        for node in meshList:
            file.write("{:21.14E} {:21.14E} {:21.14E} ".format(\
                       node[0], node[1], node[2])\
                      +"{:21.14E} {:21.14E} {:21.14E} {:2d}\n".format(\
                                 node[3], node[4], node[5], node[6]))
            continue
        pass
    return meshList

# Check if verticies are equal
def isVertexEqual(v1, v2):
    if round(v1.x, 10)==round(v2.x, 10) \
    and round(v1.y, 10)==round(v2.y, 10) \
    and round(v1.z, 10)==round(v2.z, 10):
        return True
    else:
        return False
    pass

# Class Segment child of class Line
class Segment(Line):
    def __init__(self, v1_in=Vertex(), v2_in=Vertex(), port_in=0):
        super().__init__(v1_in, v2_in, port_in)
        self.adjacents = []
        self.basis = []
        self.In = 0
        self.Out = 0
        self.index = 0
        pass

# Class Basis made of segments
class Basis:
    def __init__(self, rn_m_in=Vertex(), rn_in=Vertex(), rn_p_in=Vertex(), port_in=0):
        self.rn_m = rn_m_in
        self.rn = rn_in
        self.rn_p = rn_p_in
        self.port = port_in
        pass

# Create basis from Mesh
def createBasis(myMesh, basisFileName):
    Segments = []
    count = 0
    # Convert lines to segments
    for entry in myMesh:
        newV1 = Vertex(entry[0], entry[1], entry[2])
        newV2 = Vertex(entry[3], entry[4], entry[5])
        newSegment =  Segment(newV1, newV2, entry[6])
        count+=1
        newSegment.index = count
        Segments.append(newSegment)
        continue
    # Find segments adjecents
    for segS in Segments:
        for segD in Segments:
            if (isVertexEqual(segS.v1, segD.v1) or \
            isVertexEqual(segS.v1, segD.v2) or \
            isVertexEqual(segS.v2, segD.v1) or \
            isVertexEqual(segS.v2, segD.v2)) and \
            segS.index not in segD.adjacents and \
            segD.index not in segS.adjacents and \
            segS.index!=segD.index:
                segS.adjacents.append(segD.index)
                segD.adjacents.append(segS.index)
                pass
            continue
        continue
    logFileName = "Mesh/MeshLog.txt"
    logFile = open(logFileName, "w")
    # Logging segments adjacents
    for segment in Segments:
        for i in range(len(segment.adjacents)):
            logFile.write("Segment {} is adjacent to {}\n".format(segment.index, segment.adjacents[i]))
            continue
        continue

    BasisList = []
    maxBasis = 2*len(Segments)
    count = 0
    # Iterate over all segments
    for i in range(len(Segments)):
        segS = Segments[i]
        logFile.write("{}->".format(segS.index))
        for i in range(maxBasis):
            flag = True
            if len(segS.adjacents)>len(segS.basis):
                count+=1
                for segD_index in segS.adjacents:
                    if segD_index not in segS.basis and flag and \
                    segS.index not in Segments[segD_index-1].basis:
                        segD = Segments[segD_index-1]
                        segS.basis.append(segD.index)
                        segD.basis.append(segS.index)
                        logFile.write("{}->".format(segD_index))
                        # Check which scenario applies
                        new_S_v1 =  Vertex(segS.v1.x, segS.v1.y, segS.v1.z)
                        new_S_v2 =  Vertex(segS.v2.x, segS.v2.y, segS.v2.z)
                        new_D_v1 =  Vertex(segD.v1.x, segD.v1.y, segD.v1.z)
                        new_D_v2 =  Vertex(segD.v2.x, segD.v2.y, segD.v2.z)
                        if  isVertexEqual(new_S_v2, new_D_v1):
                            rn_m =  Vertex(new_S_v1.x, new_S_v1.y, new_S_v1.z)
                            rn =  Vertex(new_S_v2.x, new_S_v2.y, new_S_v2.z)
                            rn_p =  Vertex(new_D_v2.x, new_D_v2.y, new_D_v2.z)
                            # print("Scenario I")
                            pass
                        if  isVertexEqual(new_S_v2, new_D_v2):
                            rn_m =  Vertex(new_S_v1.x, new_S_v1.y, new_S_v1.z)
                            rn =  Vertex(new_S_v2.x, new_S_v2.y, new_S_v2.z)
                            rn_p =  Vertex(new_D_v1.x, new_D_v1.y, new_D_v1.z)
                            # print("Scenario II")
                            pass
                        if  isVertexEqual(new_S_v1, new_D_v1):
                            rn_m =  Vertex(new_S_v2.x, new_S_v2.y, new_S_v2.z)
                            rn =  Vertex(new_S_v1.x, new_S_v1.y, new_S_v1.z)
                            rn_p =  Vertex(new_D_v2.x, new_D_v2.y, new_D_v2.z)
                            # print("Scenario III")
                            pass
                        if  isVertexEqual(new_S_v1, new_D_v2):
                            rn_m =  Vertex(new_S_v2.x, new_S_v2.y, new_S_v2.z)
                            rn =  Vertex(new_S_v1.x, new_S_v1.y, new_S_v1.z)
                            rn_p =  Vertex(new_D_v1.x, new_D_v1.y, new_D_v1.z)
                            # print("Scenario IV")
                            pass
                        # Include ports
                        if (segS.port==segD.port) and segS.port!=0:
                            newPort = segS.port
                            pass
                        else:
                            newPort = 0
                            pass
                        newBasis =  Basis(rn_m, rn, rn_p, newPort)
                        BasisList.append(newBasis)
                        # Swapping
                        segS = segD
                        flag = False
                        pass
                    continue
                pass
            else:
                break
            continue
        logFile.write("\n")
        continue
    logFile.write("Total of {} iterations were used\n".format(count))
    logFile.close()
    # Write to basis file
    with open(basisFileName, "w") as file:
        for newBasis in BasisList:
            file.write("{:21.14E} {:21.14E} {:21.14E} ".format(newBasis.rn_m.x, newBasis.rn_m.y, newBasis.rn_m.z)+\
            "{:21.14E} {:21.14E} {:21.14E} ".format(newBasis.rn.x, newBasis.rn.y, newBasis.rn.z)+\
            "{:21.14E} {:21.14E} {:21.14E} ".format(newBasis.rn_p.x, newBasis.rn_p.y, newBasis.rn_p.z)+\
            "{:2d}\n".format(newBasis.port))
            continue
        pass
    # Report the number of basis
    # print("{} basis were found!".format(len(BasisList)))
    return len(BasisList)

# Create a mesh from shape
def meshShape(myShape, fileName, option):
    meshList = []
    for line in myShape.Lines:
        L = line.len()
        delta = line.delta
        assert(delta>0)
        if option==True:
            if delta>L:
                delta = L
                pass
            n = L/delta
            n = math.ceil(n)
            if n%2 != 0:
                n+=1
                pass
            delta = L/n
            for i in range(0, n):
                alpha_i = (i+0)*delta/L
                node1 = Vertex()
                node1 = line.v1+(line.v2-line.v1).scale(alpha_i)
                alpha_i = (i+1)*delta/L
                node2 = Vertex()
                node2 = line.v1+(line.v2-line.v1).scale(alpha_i)
                if line.port!=0 and (i==(n/2-1) or i==(n/2)):
                    newEntrie = [ node1.x, node1.y, node1.z, \
                                  node2.x, node2.y, node2.z, line.port ]
                    pass
                else:
                    newEntrie = [ node1.x, node1.y, node1.z, \
                                  node2.x, node2.y, node2.z, 0 ]
                    pass
                meshList.append(newEntrie)
                continue
        else:
            node1 = Vertex()
            node1 = line.v1
            node2 = Vertex()
            node2 = line.v2
            if line.port!=0:
                newEntrie = [ node1.x, node1.y, node1.z, \
                              node2.x, node2.y, node2.z, line.port ]
                pass
            else:
                newEntrie = [ node1.x, node1.y, node1.z, \
                              node2.x, node2.y, node2.z, 0 ]
                pass
            meshList.append(newEntrie)
            pass
        continue
    n = len(meshList)
    # print("{} elements were found!".format(n))
    with open(fileName, "w") as file:
        for node in meshList:
            file.write("{:21.14E} {:21.14E} {:21.14E} ".format(\
                       node[0], node[1], node[2])\
                      +"{:21.14E} {:21.14E} {:21.14E} {:2d}\n".format(\
                                 node[3], node[4], node[5], node[6]))
            continue
        pass
    return meshList

# Check if verticies are equal
def isVertexEqual(v1, v2):
    if round(v1.x, 10)==round(v2.x, 10) \
    and round(v1.y, 10)==round(v2.y, 10) \
    and round(v1.z, 10)==round(v2.z, 10):
        return True
    else:
        return False
    pass

# Class Segment child of class Line
class Segment(Line):
    def __init__(self, v1_in=Vertex(), v2_in=Vertex(), port_in=0):
        super().__init__(v1_in, v2_in, port_in)
        self.adjacents = []
        self.basis = []
        self.In = 0
        self.Out = 0
        self.index = 0
        pass

# Class Basis made of segments
class Basis:
    def __init__(self, rn_m_in=Vertex(), rn_in=Vertex(), rn_p_in=Vertex(), port_in=0):
        self.rn_m = rn_m_in
        self.rn = rn_in
        self.rn_p = rn_p_in
        self.port = port_in
        pass
