import numpy as np
import time
import csv

# Library (Optional)

# Ideal Transmission Line
class TransmissionLine():
    def __init__(self, Z0=50.0, L=1.0, vp=3.0E8):
        self.Z0 = Z0
        self.L = L
        self.vp = vp
        pass
    def getSMatrix(self, Z0, freq):
        Gamma_s = (self.Z0-Z0)/(self.Z0+Z0)
        Gamma_l = (Z0-self.Z0)/(self.Z0+Z0)
        beta = 2.0*np.pi*freq/self.vp
        Gamma_0 = Gamma_l*np.exp(-1j*2.0*beta*self.L)
        S11 = (Gamma_s+Gamma_0)/(1.0+Gamma_s*Gamma_0)
        S21 = (1.0+Gamma_l)/(1.0+Gamma_s*Gamma_0)
        S21*=(2.0*self.Z0/(self.Z0+Z0))*np.exp(-1j*beta*self.L)
        return np.matrix([[S11, S21], [S21, S11]])

# Ideal Series Impedance
class SeriesImpedance():
    def __init__(self):
        pass
    def getSMatrix(self, Z0, Z):
        S11 = Z/(Z+2.0*Z0)
        S21 = 2.0*Z0/(Z+2.0*Z0)
        return np.matrix([[S11, S21], [S21, S11]])

# Ideal Shunt Impedance
class ShuntImpedance():
    def __init__(self):
        pass
    def getSMatrix(self, Z0, Y):
        S11 = -Z0*Y/(2.0+Z0*Y)
        S21 = 2.0/(2.0+Z0*Y)
        return np.matrix([[S11, S21], [S21, S11]])

# Ideal Termination Impedance
class TerminationImpedance():
    def __init__(self):
        pass
    def getSMatrix(self, Z0, Z):
        S11 = (Z-Z0)/(Z+Z0)
        return np.matrix([[S11]])

# Ideal Open Circuit
class OpenCircuit():
    def __init__(self):
        pass
    def getSMatrix(self):
        S11 = +1.0
        return np.matrix([[S11]])

# Ideal Short Circuit
class ShortCircuit():
    def __init__(self):
        pass
    def getSMatrix(self):
        S11 = -1.0
        return np.matrix([[S11]])

# Ideal Splice
class Splice():
    def __init__(self, N=3):
        assert(N>2)
        self.N = N
        pass
    def getSMatrix(self):
        N = self.N
        SMatrix = np.zeros([N, N], dtype=np.complex128)
        for m in range(N):
            for n in range(N):
                if m==n:
                    SMatrix[m, n] = 2.0/N-1.0
                    pass
                else:
                    SMatrix[m, n] = 2.0/N
                    pass
                continue
            continue
        return SMatrix

# Timer object
class Timer():
    def __init__(self):
        self.start = 0
        self.stop = 0
        self.elapsed = 0
        pass
    # Tic the timer
    def tic(self):
        self.start = time.time()
        pass
    # Toc the timer
    def toc(self):
        self.stop = time.time()
        self.elapsed = self.stop-self.start
        if self.elapsed<60:
            print("Elapsed time is {:0.4} seconds".format(self.elapsed))
            pass
        elif self.elapsed<3600:
            print("Elapsed time is {:0.4} minutes".format(self.elapsed/60))
            pass
        else:
            print("Elapsed time is {:0.4} hours".format(self.elapsed/3600))
            pass
        pass

# Read S-parameters file (formated)
def readSParemetersFiles(fileName):
    with open(fileName, 'r') as f:
        Data = list(csv.reader(f, delimiter=","))
        pass
    return np.array(Data)

# Check dimensions consistency
def checkDimensions(Data):
    dim = []
    freq = Data[0][:, 0]
    for device in Data:
        dim.append(device.shape)
        for i in range(len(freq)):
            try:
                assert(device[i, 0]==freq[i])
                pass
            except:
                print("Error: Inconsistent Frequency Vectors!")
                exit(1)
            continue
        continue
    Ns = dim[0][0]
    for dim_ in dim:
        try:
            assert(Ns==dim_[0])
            pass
        except:
            print("Error: Inconsistent Frequency Samples!")
            exit(1)
        continue
        continue
    return freq, Ns

# Get S-matrix from iteration data
def getSMatrix(Data):
    N = int(np.sqrt((Data.shape[0]-1)/2))
    SM = np.zeros([N, N], dtype=np.complex128)
    count = 1
    for m in range(0, N):
        for n in range(0, N):
            SM[m, n] = np.complex128(Data[count])+1j*np.complex128(Data[count+1])
            count+=2
            continue
        continue
    return SM

# Save S-matrix results
def saveSMatrix(SMatrix, freq, fileName):
    Ns = SMatrix.shape[0]
    N = SMatrix.shape[1]
    with open(fileName, 'w') as file:
        for i in range(0, Ns):
            str = "{:21.14E},".format(float(freq[i]))
            file.write(str)
            for m in range(0, N):
                for n in range(0, N):
                    str = "{:21.14E},".format(np.real(SMatrix[i, m, n]))
                    file.write(str)
                    if not(m==N-1 and n==N-1):
                        str = "{:21.14E},".format(np.imag(SMatrix[i, m, n]))
                        file.write(str)
                        pass
                    else:
                        str = "{:21.14E}".format(np.imag(SMatrix[i, m, n]))
                        file.write(str)
                        pass
                    continue
                continue
            file.write("\n")
            continue
        pass
    pass

# Node object
class Node():
    def __init__(self, index=0, type='a', inputs=[]):
        self.index = index
        self.type = type
        self.inputs = inputs
        self.port = 0
        self.device = 0
        self.port_index_gbl = 0
        pass
    def __repr__(self):
        str = "Node {}{} in Port {} from Device {}".format(\
        self.type, self.index, self.port, self.device)
        str+=" is connected to: "
        for input in self.inputs:
            str+="{}{} from Device {} ({}), ".format(input[0].type,\
            input[0].index, input[0].device, input[1])
            continue
        str+="\n"
        return str

# Port object
class Port():
    def __init__(self, index=0, device=0, nodes=[]):
        self.index = index
        self.device = device
        self.nodes = nodes
        self.isConnected = False
        self.isExcited = False
        self.index_gbl = 0
        pass
    def __repr__(self):
        str = ""
        for node in self.nodes:
            str+=node.__repr__()
            continue
        return str

# Device object
class Device():
    def __init__(self, index=0, n_ports=1):
        self.index = index
        self.n_ports = n_ports
        self.SMatrix = []
        self.ports = []
        self.flag = False
        pass
    def setSMatrix(self, SMatrix):
        try:
            dim = SMatrix.shape
            assert(dim[0]==dim[1] and dim[0]==self.n_ports)
            pass
        except:
            print("Error: Wrong S-Matrix Size!")
            exit(1)
        self.SMatrix = SMatrix
        self.flag = True
        self.getPorts()
        pass
    def getPorts(self):
        self.ports = []
        for i in range(0, self.n_ports):
            newNode_a = Node(i+1, type='a')
            newNode_a.inputs = []
            newNode_a.port = i+1
            newNode_a.port = i+1
            newNode_a.device = self.index
            newNode_b = Node(i+1, type='b')
            newNode_b.inputs = []
            newNode_b.port = i+1
            newNode_b.port = i+1
            newNode_b.device = self.index
            newPort = Port(i+1, self.index)
            newPort.nodes = []
            newPort.nodes.append(newNode_a)
            newPort.nodes.append(newNode_b)
            self.ports.append(newPort)
            continue
        for portS in self.ports:
            for portD in self.ports:
                for nodeS in portS.nodes:
                    if nodeS.type=='a':
                        for nodeD in portD.nodes:
                            if nodeD.type=='b':
                                nodeD.inputs.append([nodeS, self.SMatrix[nodeD.index-1, nodeS.index-1]])
                                pass
                            continue
                        pass
                    continue
                continue
            continue
        pass
    def __repr__(self):
        str = "Device {}:\n".format(self.index)
        for port in self.ports:
            str+=port.__repr__()
            continue
        return str

# Network object
class Network():
    def __init__(self, index=0, devices=[]):
        self.index = index
        self.devices = devices
        self.SMatrix = []
        pass
    def connect(self, deviceA, portA, deviceB, portB):
        deviceA = deviceA.index
        deviceB = deviceB.index
        for device in self.devices:
            if device.index==deviceA:
                try:
                    assert(device.flag==True)
                    pass
                except:
                    print("Error: Missing S-Matrix!")
                    exit(1)
                deviceS = device
                pass
            if device.index==deviceB:
                try:
                    assert(device.flag==True)
                    pass
                except:
                    print("Error: Missing S-Matrix!")
                    exit(1)
                deviceD = device
            continue
        for port in deviceS.ports:
            if port.index==portA:
                portS = port
            continue
        for port in deviceD.ports:
            if port.index==portB:
                portD = port
            continue
        for node in portS.nodes:
            if node.type=='a':
                nodeSa = node
            if node.type=='b':
                nodeSb = node
            continue
        for node in portD.nodes:
            if node.type=='a':
                nodeDa = node
            if node.type=='b':
                nodeDb = node
            continue
        nodeDa.inputs.append([nodeSb, 1.0])
        nodeSa.inputs.append([nodeDb, 1.0])
        try:
            assert(portS!=portD)
            assert(not portS.isConnected)
            assert(not portD.isConnected)
            pass
        except:
            print("Error: Overwriting Port Connection!")
            exit(1)
        portS.isConnected = True
        portD.isConnected = True
        pass
    def __repr__(self):
        str = "Network {}: consist of {} Devices\n\n".format(self.index, len(self.devices))
        for device in self.devices:
            str+=device.__repr__()
            str+="\n"
            continue
        for device in self.devices:
            for port in device.ports:
                if port.isExcited==True:
                    str+="Port {} from Device {} is Excitation Port {} (globally)\n".format(port.index, device.index, port.index_gbl)
                    pass
                continue
            continue
        str+="\n"
        for device in self.devices:
            for port in device.ports:
                if port.isConnected==True:
                    for node in port.nodes:
                        for input in node.inputs:
                            if input[0].device!=device.index:
                                str+="Port {} from Device {} is connected to Port {} from Device {}\n".format(port.index, device.index, input[0].port, input[0].device)
                                pass
                            continue
                        continue
                    pass
                continue
            continue
        return str
    def excitation(self, deviceA, portA, index):
        for device in self.devices:
            if device==deviceA:
                for port in device.ports:
                    if port.index == portA:
                        port.index_gbl = index
                        port.isExcited = True
                        for node in port.nodes:
                            node.port_index_gbl = port.index_gbl
                            continue
                        break
                    continue
                pass
            continue
        pass
    def check(self, logFlag):
        list = []
        for device in self.devices:
            try:
                assert(device.flag==True)
                pass
            except:
                print("Error: Missing S-Matrix!")
                exit(1)
            list.append(device.index)
            continue
        try:
            assert(len(set(list))==len(list))
            pass
        except:
            print("Error: Incorrect Devices Indices!")
            exit(1)
        list = []
        countAll = 0
        for device in self.devices:
            for port in device.ports:
                countAll+=1
                if port.isExcited==True:
                    list.append(port.index_gbl)
                    pass
                if port.isExcited and port.isConnected:
                    print("Error: Duplicate Port Connection and Excitation!")
                    exit(1)
                continue
            continue
        try:
            assert(len(list)>0)
            assert(len(set(list))==len(list))
            pass
        except:
            print("Error: Incorrect Global Ports Indices!")
            exit(1)
        count = 0
        for device in self.devices:
            for port in device.ports:
                if port.isConnected==True:
                    count+=1
                    pass
                continue
            continue
        try:
            assert(count+len(list)==countAll)
            pass
        except:
            print("Error: Incorrect Ports Connections Indices!")
            exit(1)
        if logFlag:
            with open("log.txt", 'w') as file:
                file.write(self.__repr__())
                pass
            pass
        pass
    def solve(self, logFlag=False):
        self.check(logFlag)
        nodes_a = []
        nodes_b = []
        nodes_x = []
        for device in self.devices:
            for port in device.ports:
                if port.isExcited==True:
                    for node in port.nodes:
                        if node.type=='a':
                            nodes_a.append(node)
                            pass
                        if node.type=='b':
                            nodes_b.append(node)
                            pass
                        continue
                    pass
                else:
                    for node in port.nodes:
                        nodes_x.append(node)
                        continue
                    pass
                    pass
                continue
            continue
        nodes_a.sort(key=lambda x: x.port_index_gbl, reverse=False)
        nodes_b.sort(key=lambda x: x.port_index_gbl, reverse=False)
        Wba = np.zeros([len(nodes_b), len(nodes_a)], dtype=np.complex128)
        Wbx = np.zeros([len(nodes_b), len(nodes_x)], dtype=np.complex128)
        Wxa = np.zeros([len(nodes_x), len(nodes_a)], dtype=np.complex128)
        Wxx = np.zeros([len(nodes_x), len(nodes_x)], dtype=np.complex128)
        # Fill Wba and Wbx
        m = 0
        for node_b in nodes_b:
            n = 0
            for node_a in nodes_a:
                for nodeS in node_b.inputs:
                    if nodeS[0]==node_a:
                        Wba[m, n] =nodeS[1]
                        pass
                    continue
                n+=1
                continue
            n = 0
            for node_x in nodes_x:
                for nodeS in node_b.inputs:
                    if nodeS[0]==node_x:
                        Wbx[m, n] =nodeS[1]
                        pass
                    continue
                n+=1
                continue
            m+=1
            continue
        # Fill Wxa and Wxx
        m = 0
        for node_x in nodes_x:
            n = 0
            for node_a in nodes_a:
                for nodeS in node_x.inputs:
                    if nodeS[0]==node_a:
                        Wxa[m, n] =nodeS[1]
                        pass
                    continue
                n+=1
                continue
            n = 0
            for node_x_ in nodes_x:
                for nodeS in node_x.inputs:
                    if nodeS[0]==node_x_:
                        Wxx[m, n] =nodeS[1]
                        pass
                    continue
                n+=1
                continue
            m+=1
            continue
        I = np.eye(len(nodes_x))
        self.SMatrix = Wba+np.matmul(Wbx, np.matmul(np.linalg.inv(I-Wxx), Wxa))
        pass




