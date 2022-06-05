#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
//
#include "../include/TestBench.h"
#include "../include/Vector.h"
#include "../include/Utilities.h"
#include "../include/Shape.h"
#include "../include/Matrix.h"
#include "../include/QuadL.h"
#include "../include/MM_Engine.h"

void TestVector(){

  Vector A={1.0, -0.5, 1.7};
  RealVector B=DefaultRealVector;
  Vector C=subVector(A, RealVector2Vector(B));
  showVector(C);

}

void TestUtilities(){
  Timer T=DefaultTimer;
  ticTimer(&T);
  int N=100;
  LogFile myLogFile=DefaultLogFile;
  openLogFile(&myLogFile, "Data/LogFile.txt");
  char newLine[1000];
  for (int i=0; i<N; i++){
    sprintf(newLine, "Step %3d/%3d\n", i+1, N);
    appendLogFile(myLogFile, newLine);
    usleep(1000);
    progressBar(i, N);
  }
  closeLogFile(&myLogFile);
  tocTimer(&T);
  showTimer(T);

}

void TestShape(){
  double L=0.5;
  // double radius=0.5;
  double delta=0.1;
  // Test Shapes
  createDipole_z(L, delta); // Vertrical Dipole
  // createCircle_xy(radius, delta); // Horizontal Circle
  //
  Shape myShape=createShape();
  logShape(myShape, "Data/ShapeLog.txt");
  deleteShape(&myShape);
}

void TestVerticalDipole(){
  double L=0.47;
  double a=5.0E-3;
  double segPerLambda=65.0;
  // Frequency
  double lambda=1.0;
  // Create shape
  createDipole_z(L, 1.0/(segPerLambda*lambda));
  Shape myShape=createShape();
  myShape.a = a;
  printf("%d basis functions\n", myShape.NBasis);
  logShape(myShape, "Data/ShapeLog.txt");
  // Add ports
  PortsConfig myPorts=DefaultPortsConfig;
  myPorts.NExcitations = 1;
  myPorts.ExcitationPorts = malloc(myPorts.NExcitations*sizeof(int));
  myPorts.Excitations = malloc(myPorts.NExcitations*sizeof(complex double));
  myPorts.ExcitationPorts[0] = 1;
  myPorts.Excitations[0] = +1.0;
  //
  Matrix Zmn=DefaultMatrix;
  Matrix Vm=DefaultMatrix;
  Matrix In=DefaultMatrix;
  int NBasis=myShape.NBasis;
  allocateMatrix(&Zmn, NBasis, NBasis);
  allocateMatrix(&Vm, NBasis, 1);
  allocateMatrix(&In, NBasis, 1);
  //
  setQuadL(16);
  MM_DeltaGap(&Zmn, &Vm, &In, myShape, lambda, myPorts);
  unsetQuadL();
  saveMatrix(Zmn, "Data/Zmn.dat");
  GEPP(Zmn, Vm, &In);
  saveMatrix(In, "Data/In.dat");
  complex double Zin=0.0;
  int flag=0;
  for (int m=0; m<NBasis; m++){
    if (myShape.BasisList[m].port==1){
      Zin=1.0/In.data[m][0];
      flag++;
    }
  }
  assert(flag==1);
  printf("(%0.2f, %0.2f)\n", creal(Zin), cimag(Zin));
  //
  free(myPorts.ExcitationPorts);
  free(myPorts.Excitations);
  free(myPorts.LoadPorts);
  free(myPorts.Loads);
  deallocateMatrix(&Zmn);
  deallocateMatrix(&Vm);
  deallocateMatrix(&In);
  deleteShape(&myShape);
}