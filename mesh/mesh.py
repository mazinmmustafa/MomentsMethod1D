import numpy as np

class Vertex():
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x 
        self.y = y 
        self.z = z 
        pass
    def __add__(self, other):
        return Vertex(self.x+other.x, self.y+other.y, self.z+other.z)
    def __sub__(self, other):
        return Vertex(self.x-other.x, self.y-other.y, self.z-other.z)
    def __rmul__(self, other):
        return Vertex(other*self.x, other*self.y, other*self.z)
    def __lmul__(self, other):
        return Vertex(other*self.x, other*self.y, other*self.z)
    def __truediv__(self, other):
        return Vertex(self.x/other, self.y/other, self.z/other)
    def __repr__(self):
        return "<{:21.14E}, {:21.14E}, {:21.14E}>".format(self.x, self.y, self.z)
    def mag(self):
        return np.sqrt(self.x**2+self.y**2+self.z**2)
    def unit(self):
        mag_ = self.mag()
        return Vertex(self.x/mag_, self.y/mag_, self.z/mag_)
    
def is_Vertex_equal(v1, v2, tol):
    if (v1-v2).mag()<=tol:
        return True
    else:
        return False

class Segment():
    def __init__(self, v1=Vertex(), v2=Vertex(), is_port=False, port_number=0, index=0):
        self.v1 = v1
        self.v2 = v2
        try:
            assert(index>=0)
            pass
        except:
            print("Error: invalid segment index!")
            exit(1)
        self.index = index
        self.is_port = is_port
        self.port_number = port_number
        self.adjacents_list = []
        pass
    def __repr__(self):
        string = f"Segment {self.index}: Port {self.port_number}\n"
        string+= f"{self.v1}\n"
        string+= f"{self.v2}\n"
        return string

class Basis():
    def __init__(self, r_m=Vertex(), r_n=Vertex(), r_p=Vertex(), port_number=0):
        self.r_m = r_m
        self.r_n = r_n
        self.r_p = r_p
        self.port_number = port_number
        pass

class Line():
    def __init__(self, v1=Vertex(), v2=Vertex(), is_port=False, port_number=0):
        self.v1 = v1
        self.v2 = v2
        self.is_port = is_port
        try:
            assert(port_number>=0)
            if is_port:
                assert(port_number>0)
                pass
            else:
                assert(port_number==0)
                pass
            pass
        except:
            print("Error: invalid port number!")
            exit(1)
        pass
        if not is_port:
            self.port_number = 0
            pass
        else:
            self.port_number = port_number
            pass
        pass

class Shape():
    def __init__(self, lambda_=1.0, lc=0.1):
        self.lines_list = []
        self.segments_list = []
        self.bases_list = []
        try:
            assert(lambda_>0.0)
            pass
        except:
            print("Error: invalid value for lambda!")
            exit(1)
        try:
            assert(lc<1.0)
            pass
        except:
            print("Error: invalid lc value!")
            exit(1)
        self.lambda_ = lambda_
        self.lc = lc
        self.tol = 1E-6
        self.index = 0
        self.defined_ports = []
        pass
    def add_line(self, line):
        line.v1 = line.v1/self.lambda_
        line.v2 = line.v2/self.lambda_
        try:
            assert(not is_Vertex_equal(line.v1, line.v2, self.tol))
            pass
        except:
            print("Error: invalid line definition!")
            exit(1)
        self.lines_list.append(line)
        pass
    def add_segment(self, segment):
        segment.v1 = segment.v1/self.lambda_
        segment.v2 = segment.v2/self.lambda_
        try:
            assert(not is_Vertex_equal(segment.v1, segment.v2, self.tol))
            pass
        except:
            print("Error: invalid segment definition!")
            exit(1)
        pass
        segment.index = self.index
        self.index+=1
        self.segments_list.append(segment)
        pass
    def get_basis(self, segment1, segment2):
        flag = 0
        if is_Vertex_equal(segment1.v2, segment2.v1, self.tol):
            r_m = segment1.v1
            r_n = segment1.v2
            r_p = segment2.v2
            flag+=1
            pass
        if is_Vertex_equal(segment1.v2, segment2.v2, self.tol):
            r_m = segment1.v1
            r_n = segment1.v2
            r_p = segment2.v1
            flag+=1
            pass
        if is_Vertex_equal(segment1.v1, segment2.v1, self.tol):
            r_m = segment1.v2
            r_n = segment1.v1
            r_p = segment2.v2
            flag+=1
            pass
        if is_Vertex_equal(segment1.v1, segment2.v2, self.tol):
            r_m = segment1.v2
            r_n = segment1.v1
            r_p = segment2.v1
            flag+=1
            pass
        try:
            assert(flag==0 or flag==1)
            pass
        except:
            print("Error: invalid basis definition!")
            exit(1)
        pass
        if flag==1:
            if segment1.port_number==segment2.port_number and \
                segment1.is_port and segment2.is_port:
                port_number = segment1.port_number
                basis = Basis(r_m, r_n, r_p, port_number)
                if port_number>0:
                    try:
                        assert(port_number not in self.defined_ports)
                        self.defined_ports.append(port_number)
                        pass
                    except:
                        print("Error: port number was already defined!")
                        exit(1)
                    pass
                pass
            else:
                basis = Basis(r_m, r_n, r_p, 0)
                pass
            self.bases_list.append(basis)
            segment1.adjacents_list.append(segment2.index)
            segment2.adjacents_list.append(segment1.index)
            pass
        pass
    def mesh(self):
        for line in self.lines_list:
            n = int(np.ceil((line.v1-line.v2).mag()/self.lc))
            if n%2!=0:
                n+=1
                pass 
            dl = (line.v1-line.v2).mag()/n
            l_hat = (line.v2-line.v1).unit()
            for i in range(n):
                v1 = line.v1+(i+0)*dl*l_hat
                v2 = line.v1+(i+1)*dl*l_hat
                if i==n/2 or i==n/2-1:
                    segment = Segment(v1, v2, True, line.port_number, self.index)
                    pass
                else:
                    segment = Segment(v1, v2, False, 0, self.index)
                    pass
                self.segments_list.append(segment)
                self.index+=1
                continue
            continue
        
        for segment_s in self.segments_list:
            for segment_d in self.segments_list:
                if segment_s.index!=segment_d.index and \
                    segment_s.index not in segment_d.adjacents_list and \
                    segment_d.index not in segment_s.adjacents_list:
                    try:
                        assert(segment_s.index != segment_d.index)
                        pass
                    except:
                        print("Error: invalid segments indices!")
                        exit(1)
                    self.get_basis(segment_s, segment_d)
                    pass
                continue
            continue
        self.write()
        pass
    def write(self):
        file = open("mesh/shape_info.txt", "w")
        file.write("{:d}".format(len(self.bases_list)))
        file.close()
        file = open("mesh/bases_info.txt", "w")
        for basis in self.bases_list:
            file.write("{:21.14E} ".format(basis.r_m.x))
            file.write("{:21.14E} ".format(basis.r_m.y))
            file.write("{:21.14E} ".format(basis.r_m.z))
            file.write("{:21.14E} ".format(basis.r_n.x))
            file.write("{:21.14E} ".format(basis.r_n.y))
            file.write("{:21.14E} ".format(basis.r_n.z))
            file.write("{:21.14E} ".format(basis.r_p.x))
            file.write("{:21.14E} ".format(basis.r_p.y))
            file.write("{:21.14E} ".format(basis.r_p.z))
            file.write("{:2d}\n".format(basis.port_number))
            continue
        file.close()


