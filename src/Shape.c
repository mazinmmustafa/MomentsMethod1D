#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <windows.h>
//
#include "../include/Shape.h"
#include "../include/Utilities.h"

Shape createShape(){
  system("Mesh\\shape.py");
  Shape myShape=DefaultShape;
  FILE *file=NULL;
  // Read number of basis
  file = fopen("Mesh/BasisInfo.txt", "r");
  assert(file!=NULL);
  fscanf(file, "%d", &myShape.NBasis);
  assert(myShape.NBasis>0);
  fclose(file);
  // Read basis data
  myShape.BasisList = malloc(myShape.NBasis*sizeof(Basis));
  file = fopen("Mesh/Basis.dat", "r");
  assert(file!=NULL);
  double x, y, z;
  int port;
  Basis newBasis=DefaultBasis;
  for (int i=0; i<myShape.NBasis; i++){
    fscanf(file, "%lE", &x); fscanf(file, "%lE", &y); fscanf(file, "%lE", &z);
    newBasis.rn_m.x = x; newBasis.rn_m.y = y; newBasis.rn_m.z = z;
    fscanf(file, "%lE", &x); fscanf(file, "%lE", &y); fscanf(file, "%lE", &z);
    newBasis.rn.x = x; newBasis.rn.y = y; newBasis.rn.z = z;
    fscanf(file, "%lE", &x); fscanf(file, "%lE", &y); fscanf(file, "%lE", &z);
    newBasis.rn_p.x = x; newBasis.rn_p.y = y; newBasis.rn_p.z = z;
    fscanf(file, "%d", &port);
    newBasis.port = port;
    myShape.BasisList[i] = newBasis;
  }
  fclose(file);
  // Process the Basis functions
  for (int i=0; i<myShape.NBasis; i++){
    myShape.BasisList[i].Ln_m = subVector(myShape.BasisList[i].rn, myShape.BasisList[i].rn_m);
    myShape.BasisList[i].Ln_p = subVector(myShape.BasisList[i].rn_p, myShape.BasisList[i].rn);
    myShape.BasisList[i].ln_m = magVector(myShape.BasisList[i].Ln_m);
    myShape.BasisList[i].ln_p = magVector(myShape.BasisList[i].Ln_p);
  }
  // Retrun the final Shape
  return myShape;
}

void deleteShape(Shape *myShape){
  free((*myShape).BasisList);
  (*myShape).BasisList = NULL;
  (*myShape).NBasis = 0;
}

void logShape(const Shape myShape, char *fileName){
  assert(myShape.BasisList!=NULL);
  assert(myShape.NBasis>0);
  LogFile myLogFile=DefaultLogFile;
  openLogFile(&myLogFile, fileName);
  char newLine[1000];
  sprintf(newLine, "Number of basis is %d\n", myShape.NBasis);
  appendLogFile(myLogFile, newLine);
  for (int i=0; i<myShape.NBasis; i++){
    sprintf(newLine, "%2d ", myShape.BasisList[i].port);
    appendLogFile(myLogFile, newLine);
    sprintf(newLine, "%21.14E %21.14E %21.14E ",
              creal(myShape.BasisList[i].rn_m.x),
              creal(myShape.BasisList[i].rn_m.y),
              creal(myShape.BasisList[i].rn_m.z));
    appendLogFile(myLogFile, newLine);
    sprintf(newLine, "%21.14E %21.14E %21.14E ",
              creal(myShape.BasisList[i].rn.x),
              creal(myShape.BasisList[i].rn.y),
              creal(myShape.BasisList[i].rn.z));
    appendLogFile(myLogFile, newLine);
    sprintf(newLine, "%21.14E %21.14E %21.14E ",
              creal(myShape.BasisList[i].rn_p.x),
              creal(myShape.BasisList[i].rn_p.y),
              creal(myShape.BasisList[i].rn_p.z));
    appendLogFile(myLogFile, newLine);
    sprintf(newLine, "%21.14E %21.14E %21.14E ",
              creal(myShape.BasisList[i].Ln_m.x),
              creal(myShape.BasisList[i].Ln_m.y),
              creal(myShape.BasisList[i].Ln_m.z));
    appendLogFile(myLogFile, newLine);
    sprintf(newLine, "%21.14E %21.14E %21.14E ",
              creal(myShape.BasisList[i].Ln_p.x),
              creal(myShape.BasisList[i].Ln_p.y),
              creal(myShape.BasisList[i].Ln_p.z));
    appendLogFile(myLogFile, newLine);
    sprintf(newLine, "%21.14E %21.14E\n",
              myShape.BasisList[i].ln_m,
              myShape.BasisList[i].ln_p);
    appendLogFile(myLogFile, newLine);
  }
  closeLogFile(&myLogFile);
}

void createPortsConfig(PortsConfig *myPorts, int NExcitations, int NLoads){
  #define myPorts (*myPorts)
  myPorts.NExcitations = NExcitations;
  myPorts.NLoads = NLoads;
  myPorts.ExcitationPorts = malloc(myPorts.NExcitations*sizeof(int));
  myPorts.Excitations = malloc(myPorts.NExcitations*sizeof(complex double));
  myPorts.LoadPorts = malloc(myPorts.NLoads*sizeof(int));
  myPorts.Loads = malloc(myPorts.NLoads*sizeof(complex double));
  #undef myPorts
}

void deletePortsConfig(PortsConfig *myPorts){
  #define myPorts (*myPorts)
  free(myPorts.ExcitationPorts);
  free(myPorts.Excitations);
  free(myPorts.LoadPorts);
  free(myPorts.Loads);
  #undef myPorts
}

void createDipole_z(double L, double delta){
  // Open the file
  FILE *file=fopen("Mesh/shape.py", "w");
  assert(file!=NULL);
  // Start writing python script
  fprintf(file, "import Mesh1D as msh\n");
  fprintf(file, "# Definitions\n");
  fprintf(file, "basisFileName = 'Mesh/Basis.dat'\n");
  fprintf(file, "msehFileName = 'Mesh/Mesh.dat'\n");
  fprintf(file, "delta = %21.14E\n", delta);
  fprintf(file, "# Define vertices\n");
  fprintf(file, "v1 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, -L/2.0);
  fprintf(file, "v2 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, +L/2.0);
  fprintf(file, "# Define lines\n");
  fprintf(file, "Line1 = msh.Line(v1, v2, 1, delta)\n");
  fprintf(file, "# Define shape\n");
  fprintf(file, "Shape1 = msh.Shape([Line1])\n");
  fprintf(file, "# Create mesh\n");
  fprintf(file, "myMesh = msh.meshShape(Shape1, msehFileName, True)\n");
  fprintf(file, "# Create Basis\n");
  fprintf(file, "N_Basis = msh.createBasis(myMesh, basisFileName)\n");
  fprintf(file, "# Save No. of basis\n");
  fprintf(file, "with open('Mesh/BasisInfo.txt', 'w') as file:\n");
	fprintf(file, "\tfile.write('{}'.format(N_Basis))\n");
	fprintf(file, "\tpass\n");
  // Close the file
  fclose(file);
}

void createCircle_xy(double radius, double delta){
  // Open the file
  FILE *file=fopen("Mesh/shape.py", "w");
  assert(file!=NULL);
  // Start writing python script
  fprintf(file, "import Mesh1D as msh\n");
  fprintf(file, "import math\n");
  fprintf(file, "# Definitions\n");
  fprintf(file, "basisFileName = 'Mesh/Basis.dat'\n");
  fprintf(file, "msehFileName = 'Mesh/Mesh.dat'\n");
  fprintf(file, "delta = %21.14E\n", delta);
  fprintf(file, "radius = %21.14E\n", radius);
  fprintf(file, "s = 2.0*math.pi*radius\n");
  fprintf(file, "N = math.ceil(s/delta)\n");
  fprintf(file, "if N%%2 != 0:\n");
	fprintf(file, "\tN+=1\n");
	fprintf(file, "\tpass\n");
  fprintf(file, "dphi = 2.0*math.pi/N\n");
  fprintf(file, "# Define vertices and line\n");
  fprintf(file, "LinesList = []\n");
  fprintf(file, "for i in range(N):\n");
	fprintf(file,"\tphi = (i+0)*dphi\n");
	fprintf(file, "\tv1 = msh.Vertex(radius*math.cos(phi), radius*math.sin(phi), 0.0)\n");
	fprintf(file, "\tphi = (i+1)*dphi\n");
	fprintf(file, "\tv2 = msh.Vertex(radius*math.cos(phi), radius*math.sin(phi), 0.0)\n");
	fprintf(file, "\tif i==0 or i==(N-1):\n");
	fprintf(file, "\t\tLine = msh.Line(v1, v2, 1, delta)\n");
	fprintf(file, "\t\tpass\n");
	fprintf(file, "\telse:\n");
	fprintf(file, "\t\tLine = msh.Line(v1, v2, 0, delta)\n");
	fprintf(file, "\t\tpass\n");
	fprintf(file, "\tLinesList.append(Line)\n");
	fprintf(file, "\tcontinue\n");
  fprintf(file, "# Define shape\n");
  fprintf(file, "Shape1 = msh.Shape(LinesList)\n");
  fprintf(file, "# Create mesh\n");
  fprintf(file, "myMesh = msh.meshShape(Shape1, msehFileName, False)\n");
  fprintf(file, "# Create Basis\n");
  fprintf(file, "N_Basis = msh.createBasis(myMesh, basisFileName)\n");
  fprintf(file, "# Save No. of basis\n");
  fprintf(file, "with open('Mesh/BasisInfo.txt', 'w') as file:\n");
	fprintf(file, "\tfile.write('{}'.format(N_Basis))\n");
	fprintf(file, "\tpass\n");
  // Close the file
  fclose(file);
}

void createTL(double L, double S, double delta){
  // Open the file
  FILE *file=fopen("Mesh/shape.py", "w");
  assert(file!=NULL);
  // Start writing python script
  fprintf(file, "import Mesh1D as msh\n");
  fprintf(file, "# Definitions\n");
  fprintf(file, "basisFileName = 'Mesh/Basis.dat'\n");
  fprintf(file, "msehFileName = 'Mesh/Mesh.dat'\n");
  fprintf(file, "delta = %21.14E\n", delta);
  fprintf(file, "# Define vertices\n");
  fprintf(file, "v1 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", -L/2.0, -S/2.0, 0.0);
  fprintf(file, "v2 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", -L/2.0, +S/2.0, 0.0);
  fprintf(file, "v3 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", +L/2.0, +S/2.0, 0.0);
  fprintf(file, "v4 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", +L/2.0, -S/2.0, 0.0);
  fprintf(file, "# Define lines\n");
  fprintf(file, "Line1 = msh.Line(v1, v2, 1, %21.14E)\n", S/3.0);
  fprintf(file, "Line2 = msh.Line(v2, v3, 0, delta)\n");
  fprintf(file, "Line3 = msh.Line(v3, v4, 2, %21.14E)\n", S/3.0);
  fprintf(file, "Line4 = msh.Line(v4, v1, 0, delta)\n");
  fprintf(file, "# Define shape\n");
  fprintf(file, "Shape1 = msh.Shape([Line1, Line2, Line3, Line4])\n");
  fprintf(file, "# Create mesh\n");
  fprintf(file, "myMesh = msh.meshShape(Shape1, msehFileName, True)\n");
  fprintf(file, "# Create Basis\n");
  fprintf(file, "N_Basis = msh.createBasis(myMesh, basisFileName)\n");
  fprintf(file, "# Save No. of basis\n");
  fprintf(file, "with open('Mesh/BasisInfo.txt', 'w') as file:\n");
	fprintf(file, "\tfile.write('{}'.format(N_Basis))\n");
	fprintf(file, "\tpass\n");
  // Close the file
  fclose(file);
}

void createGullAntenna(double h1, double h2, double h3,
               double alpha, double delta){
  alpha = deg2rad(alpha);
  // Open the file
  FILE *file=fopen("Mesh/shape.py", "w");
  assert(file!=NULL);
  // Start writing python script
  fprintf(file, "import Mesh1D as msh\n");
  fprintf(file, "# Definitions\n");
  fprintf(file, "basisFileName = 'Mesh/Basis.dat'\n");
  fprintf(file, "msehFileName = 'Mesh/Mesh.dat'\n");
  fprintf(file, "delta = %21.14E\n", delta);
  fprintf(file, "# Define vertices\n");
  fprintf(file, "v1 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", h2*sin(alpha), -h1-h2*cos(alpha)-h3, 0.0);
  fprintf(file, "v2 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", h2*sin(alpha), -h1-h2*cos(alpha), 0.0);
  fprintf(file, "v3 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, -h1, 0.0);
  fprintf(file, "v4 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, +h1, 0.0);
  fprintf(file, "v5 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", h2*sin(alpha), +h1+h2*cos(alpha), 0.0);
  fprintf(file, "v6 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", h2*sin(alpha), +h1+h2*cos(alpha)+h3, 0.0);
  fprintf(file, "# Define lines\n");
  fprintf(file, "Line1 = msh.Line(v1, v2, 0, delta)\n");
  fprintf(file, "Line2 = msh.Line(v2, v3, 0, delta)\n");
  fprintf(file, "Line3 = msh.Line(v3, v4, 1, delta)\n");
  fprintf(file, "Line4 = msh.Line(v4, v5, 0, delta)\n");
  fprintf(file, "Line5 = msh.Line(v5, v6, 0, delta)\n");
  fprintf(file, "# Define shape\n");
  fprintf(file, "Shape1 = msh.Shape([Line1, Line2, Line3, Line4, Line5])\n");
  fprintf(file, "# Create mesh\n");
  fprintf(file, "myMesh = msh.meshShape(Shape1, msehFileName, True)\n");
  fprintf(file, "# Create Basis\n");
  fprintf(file, "N_Basis = msh.createBasis(myMesh, basisFileName)\n");
  fprintf(file, "# Save No. of basis\n");
  fprintf(file, "with open('Mesh/BasisInfo.txt', 'w') as file:\n");
	fprintf(file, "\tfile.write('{}'.format(N_Basis))\n");
	fprintf(file, "\tpass\n");
  // Close the file
  fclose(file);
}

void createUdaYagiAntenna(double L1, double L2, double L3,
               double S1, double S2, double S3, int N, double delta){
  // Open the file
  FILE *file=fopen("Mesh/shape.py", "w");
  assert(file!=NULL);
  // Start writing python script
  fprintf(file, "import Mesh1D as msh\n");
  fprintf(file, "# Definitions\n");
  fprintf(file, "basisFileName = 'Mesh/Basis.dat'\n");
  fprintf(file, "msehFileName = 'Mesh/Mesh.dat'\n");
  fprintf(file, "delta = %21.14E\n", delta);
  fprintf(file, "# Define vertices\n");
  fprintf(file, "v1 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, -L1/2.0);
  fprintf(file, "v2 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, +L1/2.0);
  fprintf(file, "v3 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", S1, 0.0, -L2/2.0);
  fprintf(file, "v4 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", S1, 0.0, +L2/2.0);
  int count=5;
  for (int i=0; i<N; i++){
    fprintf(file, "v%d = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", count, S1+S2+i*S3, 0.0, -L3/2.0);
    count++;
    fprintf(file, "v%d = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", count, S1+S2+i*S3, 0.0, +L3/2.0);
    count++;
  }
  fprintf(file, "# Define lines\n");
  fprintf(file, "Lines = []\n");
  count = 1;
  for (int i=0; i<(N+2); i++){
    fprintf(file, "Line%d = msh.Line(v%d, v%d, %d, delta)\n", i+1, count, count+1, i+1);
    count+=2;
    fprintf(file, "Lines.append(Line%d)\n", i+1);
  }
  fprintf(file, "# Define shape\n");
  fprintf(file, "Shape1 = msh.Shape(Lines)\n");
  fprintf(file, "# Create mesh\n");
  fprintf(file, "myMesh = msh.meshShape(Shape1, msehFileName, True)\n");
  fprintf(file, "# Create Basis\n");
  fprintf(file, "N_Basis = msh.createBasis(myMesh, basisFileName)\n");
  fprintf(file, "# Save No. of basis\n");
  fprintf(file, "with open('Mesh/BasisInfo.txt', 'w') as file:\n");
	fprintf(file, "\tfile.write('{}'.format(N_Basis))\n");
	fprintf(file, "\tpass\n");
  // Close the file
  fclose(file);
}

void createBentWire(double L, double S, double delta){
  // Open the file
  FILE *file=fopen("Mesh/shape.py", "w");
  assert(file!=NULL);
  // Start writing python script
  fprintf(file, "import Mesh1D as msh\n");
  fprintf(file, "# Definitions\n");
  fprintf(file, "basisFileName = 'Mesh/Basis.dat'\n");
  fprintf(file, "msehFileName = 'Mesh/Mesh.dat'\n");
  fprintf(file, "delta = %21.14E\n", delta);
  fprintf(file, "# Define vertices\n");
  fprintf(file, "v1 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, -S);
  fprintf(file, "v2 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, +S);
  fprintf(file, "v3 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", L, 0.0, +S);
  fprintf(file, "v4 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", L, L, +S);
  fprintf(file, "v5 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", L, L, -S);
  fprintf(file, "v6 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", L, 0.0, -S);
  fprintf(file, "# Define lines\n");
  fprintf(file, "Line1 = msh.Line(v1, v2, 1, delta)\n");
  fprintf(file, "Line2 = msh.Line(v2, v3, 0, delta)\n");
  fprintf(file, "Line3 = msh.Line(v3, v4, 0, delta)\n");
  fprintf(file, "Line4 = msh.Line(v4, v5, 0, delta)\n");
  fprintf(file, "Line5 = msh.Line(v5, v6, 0, delta)\n");
  fprintf(file, "Line6 = msh.Line(v6, v1, 0, delta)\n");
  fprintf(file, "# Define shape\n");
  fprintf(file, "Shape1 = msh.Shape([Line1, Line2, Line3, Line4, Line5, Line6])\n");
  fprintf(file, "# Create mesh\n");
  fprintf(file, "myMesh = msh.meshShape(Shape1, msehFileName, True)\n");
  fprintf(file, "# Create Basis\n");
  fprintf(file, "N_Basis = msh.createBasis(myMesh, basisFileName)\n");
  fprintf(file, "# Save No. of basis\n");
  fprintf(file, "with open('Mesh/BasisInfo.txt', 'w') as file:\n");
	fprintf(file, "\tfile.write('{}'.format(N_Basis))\n");
	fprintf(file, "\tpass\n");
  // Close the file
  fclose(file);
}

void createBentWireSplit(double L, double S, double delta){
  // Open the file
  FILE *file=fopen("Mesh/shape.py", "w");
  assert(file!=NULL);
  // Start writing python script
  fprintf(file, "import Mesh1D as msh\n");
  fprintf(file, "# Definitions\n");
  fprintf(file, "basisFileName = 'Mesh/Basis.dat'\n");
  fprintf(file, "msehFileName = 'Mesh/Mesh.dat'\n");
  fprintf(file, "delta = %21.14E\n", delta);
  fprintf(file, "# Define vertices\n");
  fprintf(file, "v1 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, -S);
  fprintf(file, "v2 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, +S);
  fprintf(file, "v3 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", L, 0.0, +S);
  fprintf(file, "v4 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", L, L, +S);
  fprintf(file, "v5 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", L, L, -S);
  fprintf(file, "v6 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", L, 0.0, -S);
  fprintf(file, "# Define lines\n");
  fprintf(file, "Line1 = msh.Line(v1, v2, 1, delta)\n");
  fprintf(file, "Line2 = msh.Line(v2, v3, 0, delta)\n");
  fprintf(file, "Line3 = msh.Line(v3, v4, 0, delta)\n");
  fprintf(file, "Line4 = msh.Line(v4, v5, 0, delta)\n");
  fprintf(file, "Line5 = msh.Line(v5, v6, 0, delta)\n");
  fprintf(file, "Line6 = msh.Line(v6, v1, 0, delta)\n");
  fprintf(file, "Line7 = msh.Line(v3, v6, 0, delta)\n");
  fprintf(file, "# Define shape\n");
  fprintf(file, "Shape1 = msh.Shape([Line1, Line2, Line3, Line4, Line5, Line6, Line7])\n");
  fprintf(file, "# Create mesh\n");
  fprintf(file, "myMesh = msh.meshShape(Shape1, msehFileName, True)\n");
  fprintf(file, "# Create Basis\n");
  fprintf(file, "N_Basis = msh.createBasis(myMesh, basisFileName)\n");
  fprintf(file, "# Save No. of basis\n");
  fprintf(file, "with open('Mesh/BasisInfo.txt', 'w') as file:\n");
	fprintf(file, "\tfile.write('{}'.format(N_Basis))\n");
	fprintf(file, "\tpass\n");
  // Close the file
  fclose(file);
}

void createInvertedF(double L, double S, double h, double delta){
  // Open the file
  FILE *file=fopen("Mesh/shape.py", "w");
  assert(file!=NULL);
  // Start writing python script
  fprintf(file, "import Mesh1D as msh\n");
  fprintf(file, "# Definitions\n");
  fprintf(file, "basisFileName = 'Mesh/Basis.dat'\n");
  fprintf(file, "msehFileName = 'Mesh/Mesh.dat'\n");
  fprintf(file, "delta = %21.14E\n", delta);
  fprintf(file, "# Define vertices\n");
  fprintf(file, "v1 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, -h);
  fprintf(file, "v2 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, +h);
  fprintf(file, "v3 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", -S, 0.0, +h);
  fprintf(file, "v4 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", -S, 0.0, -h);
  fprintf(file, "v5 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", +L, 0.0, -h);
  fprintf(file, "v6 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", +L, 0.0, +h);
  fprintf(file, "# Define lines\n");
  fprintf(file, "Line1 = msh.Line(v1, v2, 1, delta)\n");
  fprintf(file, "Line2 = msh.Line(v2, v3, 0, delta)\n");
  fprintf(file, "Line3 = msh.Line(v3, v4, 0, delta)\n");
  fprintf(file, "Line4 = msh.Line(v4, v1, 0, delta)\n");
  fprintf(file, "Line5 = msh.Line(v1, v5, 0, delta)\n");
  fprintf(file, "Line6 = msh.Line(v6, v2, 0, delta)\n");
  fprintf(file, "# Define shape\n");
  fprintf(file, "Shape1 = msh.Shape([Line1, Line2, Line3, Line4, Line5, Line6])\n");
  fprintf(file, "# Create mesh\n");
  fprintf(file, "myMesh = msh.meshShape(Shape1, msehFileName, True)\n");
  fprintf(file, "# Create Basis\n");
  fprintf(file, "N_Basis = msh.createBasis(myMesh, basisFileName)\n");
  fprintf(file, "# Save No. of basis\n");
  fprintf(file, "with open('Mesh/BasisInfo.txt', 'w') as file:\n");
	fprintf(file, "\tfile.write('{}'.format(N_Basis))\n");
	fprintf(file, "\tpass\n");
  // Close the file
  fclose(file);
}

void createBiconicalAntenna(double h1, double h2, double h3, double S,
            int N, double delta){
  assert(N>0);
  if (N%2!=0){N++;};
  // Open the file
  FILE *file=fopen("Mesh/shape.py", "w");
  assert(file!=NULL);
  // Start writing python script
  fprintf(file, "import Mesh1D as msh\n");
  fprintf(file, "import math\n");
  fprintf(file, "# Definitions\n");
  fprintf(file, "basisFileName = 'Mesh/Basis.dat'\n");
  fprintf(file, "msehFileName = 'Mesh/Mesh.dat'\n");
  fprintf(file, "delta = %21.14E\n", delta);
  fprintf(file, "# Define vertices\n");
  fprintf(file, "v1 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, +h1);
  fprintf(file, "v2 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, +h1+h2);
  fprintf(file, "v3 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, +h1+h2+h3);
  fprintf(file, "v4 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, -h1);
  fprintf(file, "v5 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, -h1-h2);
  fprintf(file, "v6 = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", 0.0, 0.0, -h1-h2-h3);
  double dtheta=2.0*M_PI/N;
  fprintf(file, "# Define lines\n");
  fprintf(file, "Lines = []\n");
  fprintf(file, "Lines.append(msh.Line(v4, v1, 1, delta/3.0))\n");
  fprintf(file, "Lines.append(msh.Line(v1, v2, 0, delta))\n");
  fprintf(file, "Lines.append(msh.Line(v2, v3, 0, delta))\n");
  fprintf(file, "Lines.append(msh.Line(v5, v4, 0, delta))\n");
  fprintf(file, "Lines.append(msh.Line(v6, v5, 0, delta))\n");
  int count=6;
  for (int n=0; n<N; n++){
    count++;
    fprintf(file, "v%d = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", count, S*cos(n*dtheta), S*sin(n*dtheta), +h1+h2);
    fprintf(file, "Lines.append(msh.Line(v1, v%d, 0, delta))\n", count);
    fprintf(file, "Lines.append(msh.Line(v%d, v3, 0, delta))\n", count);
    // if (roundn(n*dtheta, 6)==0){
    //   fprintf(file, "Lines.append(msh.Line(v2, v%d, 0, delta))\n", count);
    // }
    count++;
    fprintf(file, "v%d = msh.Vertex( %21.14E, %21.14E, %21.14E)\n", count, S*cos(n*dtheta), S*sin(n*dtheta), -h1-h2);
    fprintf(file, "Lines.append(msh.Line(v4, v%d, 0, delta))\n", count);
    fprintf(file, "Lines.append(msh.Line(v%d, v6, 0, delta))\n", count);
    // if (roundn(n*dtheta, 6)==roundn(M_PI, 6)){
  	// if (roundn(n*dtheta, 6)==0){
    //   fprintf(file, "Lines.append(msh.Line(v5, v%d, 0, delta))\n", count);
    // }
  }
  fprintf(file, "# Define shape\n");
  fprintf(file, "Shape1 = msh.Shape(Lines)\n");
  fprintf(file, "# Create mesh\n");
  fprintf(file, "myMesh = msh.meshShape(Shape1, msehFileName, True)\n");
  fprintf(file, "# Create Basis\n");
  fprintf(file, "N_Basis = msh.createBasis(myMesh, basisFileName)\n");
  fprintf(file, "# Save No. of basis\n");
  fprintf(file, "with open('Mesh/BasisInfo.txt', 'w') as file:\n");
	fprintf(file, "\tfile.write('{}'.format(N_Basis))\n");
	fprintf(file, "\tpass\n");
  // Close the file
  fclose(file);
}