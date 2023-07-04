import rfspice as rf
import numpy as np

# Read S-parameters files for devices
Device1 = rf.readSParemetersFiles("../data.csv")

# Obtain the frequnecy vector and number of iterations
freq, Ns = rf.checkDimensions([Device1])

# Declare the overall S-matrix
SMatrix = np.zeros([Ns, 8, 8], dtype=np.complex128)

# Define ideal devices
SR1 = rf.SeriesImpedance()
SR2 = rf.SeriesImpedance()
SR3 = rf.SeriesImpedance()
SR4 = rf.SeriesImpedance()
SR5 = rf.SeriesImpedance()
SR6 = rf.SeriesImpedance()
SR7 = rf.SeriesImpedance()
SR8 = rf.SeriesImpedance()

# Set timer
T = rf.Timer()
T.tic()

Z0 = 50

# Loop over the iterations
for i in range(0, Ns):

    # Get S-matrix for each iteration
    SM1 = rf.getSMatrix(Device1[i, :])

    # Define devices: Index, Number of ports
    d1 = rf.Device(1, 8)

    d2 = rf.Device(2, 2)
    d3 = rf.Device(3, 2)
    d4 = rf.Device(4, 2)
    d5 = rf.Device(5, 2)

    d6 = rf.Device(6, 2)
    d7 = rf.Device(7, 2)
    d8 = rf.Device(8, 2)
    d9 = rf.Device(9, 2)

    pF = 1.0E-12
    j = 1j
    C_tune = 25.0*pF
    omega = 2.0*np.pi*float(freq[i])
    Z_tune = 1.0/(j*omega*C_tune)

    # Set devies S-matrices
    list = []
    d1.setSMatrix(SM1)
    list.append(d1)
    d2.setSMatrix(SR1.getSMatrix(Z0, Z_tune))
    list.append(d2)
    d3.setSMatrix(SR2.getSMatrix(Z0, Z_tune))
    list.append(d3)
    d4.setSMatrix(SR3.getSMatrix(Z0, Z_tune))
    list.append(d4)
    d5.setSMatrix(SR4.getSMatrix(Z0, Z_tune))
    list.append(d5)
    d6.setSMatrix(SR5.getSMatrix(Z0, Z_tune))
    list.append(d6)
    d7.setSMatrix(SR6.getSMatrix(Z0, Z_tune))
    list.append(d7)
    d8.setSMatrix(SR7.getSMatrix(Z0, Z_tune))
    list.append(d8)
    d9.setSMatrix(SR8.getSMatrix(Z0, Z_tune))
    list.append(d9)

    # Create a network of devices
    n1 = rf.Network(1, list)

    # Add connections between ports
    n1.connect(d1, 1, d2, 2)
    n1.connect(d1, 2, d3, 2)
    n1.connect(d1, 3, d4, 2)
    n1.connect(d1, 4, d5, 2)
    n1.connect(d1, 5, d6, 2)
    n1.connect(d1, 6, d7, 2)
    n1.connect(d1, 7, d8, 2)
    n1.connect(d1, 8, d9, 2)

    # Define excitation ports: Device, Local Port, Global Port
    n1.excitation(d2, 1, 1)
    n1.excitation(d3, 1, 2)
    n1.excitation(d4, 1, 3)
    n1.excitation(d5, 1, 4)
    n1.excitation(d6, 1, 5)
    n1.excitation(d7, 1, 6)
    n1.excitation(d8, 1, 7)
    n1.excitation(d9, 1, 8)

    # Solve the system
    n1.solve()

    # Update the results S-matrix
    SMatrix[i, :, :] = n1.SMatrix

    continue

# Unser timer
T.toc()

# Save results into data files: final S-matrix, filename
rf.saveSMatrix(SMatrix, freq, "Data/ResultsData.csv")











