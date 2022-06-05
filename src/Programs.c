#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
//
#include "../include/Programs.h"
#include "../include/Vector.h"
#include "../include/Utilities.h"
#include "../include/Shape.h"
#include "../include/Matrix.h"
#include "../include/QuadL.h"
#include "../include/MM_Engine.h"
#include "../include/FarField.h"

void programDipoleDeltaGapInputAdmittance(){
  // Definitions
  int Ns=401;
  double L_min=0.1;
  double L_max=2.0;
  double aFactor=2.0*74.2;
  double lambda=1.0;
  double segPerLambda=31.0;
  setQuadL(32);
  //
  Timer T=DefaultTimer;
  ticTimer(&T);
  double L, dL;
  double a;
  dL = (L_max-L_min)/(Ns-1.0);
  FILE *file=fopen("Data/VerticalDipole/Yin.dat", "w");
  for (int i=0; i<Ns; i++){
    L = L_min+i*dL;
    a = L/aFactor;
    //
    createDipole_z(L, lambda/segPerLambda);
    Shape myShape=createShape();
    myShape.a = a;
    // Add ports
    PortsConfig myPorts=DefaultPortsConfig;
    createPortsConfig(&myPorts, 1, 0);
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
    MM_DeltaGap(&Zmn, &Vm, &In, myShape, lambda, myPorts);
    GEPP(Zmn, Vm, &In);
    complex double Zin=0.0;
    int flag=0;
    for (int m=0; m<NBasis; m++){
      if (myShape.BasisList[m].port==1){
        Zin=1.0/In.data[m][0];
        flag++;
      }
    }
    assert(flag==1);
    fprintf(file, "%21.14E %21.14E %21.14E\n", L, creal(1.0/Zin), cimag(1.0/Zin));
    //
    deletePortsConfig(&myPorts);
    deallocateMatrix(&Zmn);
    deallocateMatrix(&Vm);
    deallocateMatrix(&In);
    deleteShape(&myShape);
    //
    progressBar(i, Ns);
  }
  fclose(file);
  unsetQuadL();
  tocTimer(&T);
  showTimer(T);
}

void programCircularLoopDeltaGapInputImpedance(){
  // Constants
  double pi=M_PI;
  // Definitions
  int Ns=401;
  double S_min=0.1;
  double S_max=2.5;
  double a=1.0E-4;
  double lambda=1.0;
  double segPerLambda=31.0;
  setQuadL(32);
  //
  Timer T=DefaultTimer;
  ticTimer(&T);
  double radius, S, dS;
  dS = (S_max-S_min)/(Ns-1.0);
  FILE *file=fopen("Data/CircularLoop/Zin.dat", "w");
  for (int i=0; i<Ns; i++){
    S = S_min+i*dS;
    radius = S/(2.0*pi);
    //
    createCircle_xy(radius, lambda/segPerLambda);
    Shape myShape=createShape();
    myShape.a = a;
    // Add ports
    PortsConfig myPorts=DefaultPortsConfig;
    createPortsConfig(&myPorts, 1, 0);
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
    MM_DeltaGap(&Zmn, &Vm, &In, myShape, lambda, myPorts);
    GEPP(Zmn, Vm, &In);
    complex double Zin=0.0;
    int flag=0;
    for (int m=0; m<NBasis; m++){
      if (myShape.BasisList[m].port==1){
        Zin=1.0/In.data[m][0];
        flag++;
      }
    }
    assert(flag==1);
    fprintf(file, "%21.14E %21.14E %21.14E\n", S, creal(Zin), cimag(Zin));
    //
    deletePortsConfig(&myPorts);
    deallocateMatrix(&Zmn);
    deallocateMatrix(&Vm);
    deallocateMatrix(&In);
    deleteShape(&myShape);
    //
    progressBar(i, Ns);
  }
  fclose(file);
  unsetQuadL();
  tocTimer(&T);
  showTimer(T);
}

void programTLDeltaGapInputImpedance(){
  // Constants
  double c0=3.0E8;
  // Definitions
  int Ns=401;
  double freq_min=1.0E5;
  double freq_max=1.0E8;
  double L=2.5;
  double S=0.03;
  double a=1.0E-3;
  double segPerLambda=31.0;
  setQuadL(32);
  //
  Timer T=DefaultTimer;
  ticTimer(&T);
  double lambda;
  double freq, dfreq;
  dfreq = (log10(freq_max)-log10(freq_min))/(Ns-1.0);
  FILE *file=fopen("Data/LoadedTL/Zin.dat", "w");
  for (int i=0; i<Ns; i++){
    freq = pow(10.0, log10(freq_min)+i*dfreq);
    lambda = c0/freq;
    //
    createTL(L, S, lambda/segPerLambda);
    Shape myShape=createShape();
    myShape.a = a;
    // Add ports
    PortsConfig myPorts=DefaultPortsConfig;
    createPortsConfig(&myPorts, 1, 1);
    myPorts.ExcitationPorts[0] = 1;
    myPorts.LoadPorts[0] = 2;
    myPorts.Excitations[0] = +1.0;
    myPorts.Loads[0] = 200.0;
    //
    Matrix Zmn=DefaultMatrix;
    Matrix Vm=DefaultMatrix;
    Matrix In=DefaultMatrix;
    int NBasis=myShape.NBasis;
    allocateMatrix(&Zmn, NBasis, NBasis);
    allocateMatrix(&Vm, NBasis, 1);
    allocateMatrix(&In, NBasis, 1);
    //
    MM_DeltaGap(&Zmn, &Vm, &In, myShape, lambda, myPorts);
    GEPP(Zmn, Vm, &In);
    complex double Zin=0.0;
    int flag=0;
    for (int m=0; m<NBasis; m++){
      if (myShape.BasisList[m].port==1){
        Zin=1.0/In.data[m][0];
        flag++;
      }
    }
    assert(flag==1);
    fprintf(file, "%21.14E %21.14E %21.14E\n", freq, creal(Zin), cimag(Zin));
    //
    deletePortsConfig(&myPorts);
    deallocateMatrix(&Zmn);
    deallocateMatrix(&Vm);
    deallocateMatrix(&In);
    deleteShape(&myShape);
    //
    progressBar(i, Ns);
  }
  fclose(file);
  unsetQuadL();
  tocTimer(&T);
  showTimer(T);
}

void programGullDeltaGapRadiationPattern(){
  // Definitions
  int Ns=1001;
  double h1=0.0714;
  double h2=0.4286;
  double h3=0.25;
  double a=5.0E-3;
  double alpha=50.0;
  double lambda=1.0;
  double segPerLambda=31.0;
  setQuadL(32);
  //
  double theta_min=-180.0;
  double theta_max=+180.0;
  double phi_min=+0.0;
  double phi_max=+360.0;
  //
  Timer T=DefaultTimer;
  ticTimer(&T);
  //
  createGullAntenna(h1, h2, h3, alpha, lambda/segPerLambda);
  Shape myShape=createShape();
  myShape.a = a;
  // Add ports
  PortsConfig myPorts=DefaultPortsConfig;
  createPortsConfig(&myPorts, 1, 0);
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
  MM_DeltaGap(&Zmn, &Vm, &In, myShape, lambda, myPorts);
  GEPP(Zmn, Vm, &In);
  FILE *file=fopen("Data/GullAntenna/In.dat", "w");
  for (int n=0; n<NBasis; n++){
    fprintf(file, "%21.14E\n", cabs(In.data[n][0]));
  }
  fclose(file);
  //
  file = fopen("Data/GullAntenna/RadiationPattern.dat", "w");
  double theta, dtheta;
  double phi, dphi;
  dtheta = (theta_max-theta_min)/(Ns-1.0);
  dphi = (phi_max-phi_min)/(Ns-1.0);
  //
  double sigmaTheta, sigmaPhi;
  for (int i=0; i<Ns; i++){
    // Azimuth
    theta = +90.0;
    phi = phi_min+i*dphi;
    sigma(myShape, theta, phi, In, lambda, &sigmaTheta, &sigmaPhi);
    fprintf(file, "%21.14E %21.14E %21.14E ", theta, phi, sigmaPhi);
    // Elevation
    phi = 0.0;
    theta = theta_min+i*dtheta;
    sigma(myShape, theta, phi, In, lambda, &sigmaTheta, &sigmaPhi);
    fprintf(file, "%21.14E %21.14E %21.14E ", theta, phi, sigmaPhi);
    //
    fprintf(file, "\n");
    progressBar(i, Ns);
  }
  fclose(file);
  //
  deletePortsConfig(&myPorts);
  deallocateMatrix(&Zmn);
  deallocateMatrix(&Vm);
  deallocateMatrix(&In);
  deleteShape(&myShape);
  unsetQuadL();
  //
  tocTimer(&T);
  showTimer(T);
}

void programUdaYagiDeltaGapRadiationPattern(){
  // Definitions
  int Ns=1001;
  double L1=0.5;
  double L2=0.47;
  double L3=0.406;
  double S1=0.25;
  double S2=0.34;
  double S3=0.34;
  int NElements=13;
  double a=3.0E-3;
  double lambda=1.0;
  double segPerLambda=101.0;
  setQuadL(32);
  //
  double theta_min=-180.0;
  double theta_max=+180.0;
  double phi_min=+0.0;
  double phi_max=+360.0;
  //
  Timer T=DefaultTimer;
  ticTimer(&T);
  //
  createUdaYagiAntenna(L1, L2, L3, S1, S2, S3, NElements, lambda/segPerLambda);
  Shape myShape=createShape();
  myShape.a = a;
  // Add ports
  PortsConfig myPorts=DefaultPortsConfig;
  createPortsConfig(&myPorts, 1, 0);
  myPorts.ExcitationPorts[0] = 2;
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
  MM_DeltaGap(&Zmn, &Vm, &In, myShape, lambda, myPorts);
  GEPP(Zmn, Vm, &In);
  FILE *file=fopen("Data/UdaYagi/In.dat", "w");
  int count=1;
  for (int n=0; n<NBasis; n++){
    if (myShape.BasisList[n].port==count){
      fprintf(file, "%21.14E\n", cabs(In.data[n][0]));
      count++;
    }
  }
  fclose(file);
  //
  file = fopen("Data/UdaYagi/RadiationPattern.dat", "w");
  double theta, dtheta;
  double phi, dphi;
  dtheta = (theta_max-theta_min)/(Ns-1.0);
  dphi = (phi_max-phi_min)/(Ns-1.0);
  //
  double sigmaTheta, sigmaPhi;
  for (int i=0; i<Ns; i++){
    // Azimuth
    theta = +90.0;
    phi = phi_min+i*dphi;
    sigma(myShape, theta, phi, In, lambda, &sigmaTheta, &sigmaPhi);
    fprintf(file, "%21.14E %21.14E %21.14E ", theta, phi, sigmaTheta);
    // Elevation
    phi = 0.0;
    theta = theta_min+i*dtheta;
    sigma(myShape, theta, phi, In, lambda, &sigmaTheta, &sigmaPhi);
    fprintf(file, "%21.14E %21.14E %21.14E ", theta, phi, sigmaTheta);
    //
    fprintf(file, "\n");
    progressBar(i, Ns);
  }
  fclose(file);
  //
  deletePortsConfig(&myPorts);
  deallocateMatrix(&Zmn);
  deallocateMatrix(&Vm);
  deallocateMatrix(&In);
  deleteShape(&myShape);
  unsetQuadL();
  //
  tocTimer(&T);
  showTimer(T);
}

void programBentWireDeltaGapInputImpedance(){
  // Constants
  double c0=3.0E8;
  // Definitions
  int Ns=1601;
  double GHz=1.0E9;
  double cm=1.0E-2;
  double freq_min=0.05*GHz;
  double freq_max=2.0*GHz;
  double L=10*cm;
  double S=3*cm;
  double a=0.5*(1.6E-3);
  double segPerLambda=31.0;
  setQuadL(32);
  //
  Timer T=DefaultTimer;
  ticTimer(&T);
  double lambda;
  double freq, dfreq;
  dfreq = (log10(freq_max)-log10(freq_min))/(Ns-1.0);
  FILE *file=fopen("Data/BentWire/Zin.dat", "w");
  for (int i=0; i<Ns; i++){
    freq = pow(10.0, log10(freq_min)+i*dfreq);
    lambda = c0/freq;
    //
    // createBentWire(L, S, lambda/segPerLambda); // Scenario I
    createBentWireSplit(L, S, lambda/segPerLambda); // Scenario II
    Shape myShape=createShape();
    myShape.a = a;
    // Add ports
    PortsConfig myPorts=DefaultPortsConfig;
    createPortsConfig(&myPorts, 1, 0);
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
    MM_DeltaGap(&Zmn, &Vm, &In, myShape, lambda, myPorts);
    GEPP(Zmn, Vm, &In);
    complex double Zin=0.0;
    int flag=0;
    for (int m=0; m<NBasis; m++){
      if (myShape.BasisList[m].port==1){
        Zin=1.0/In.data[m][0];
        flag++;
      }
    }
    assert(flag==1);
    fprintf(file, "%21.14E %21.14E %21.14E\n", freq, creal(Zin), cimag(Zin));
    //
    deletePortsConfig(&myPorts);
    deallocateMatrix(&Zmn);
    deallocateMatrix(&Vm);
    deallocateMatrix(&In);
    deleteShape(&myShape);
    //
    progressBar(i, Ns);
  }
  fclose(file);
  unsetQuadL();
  tocTimer(&T);
  showTimer(T);
}

void programBentWireDeltaGapRadiationPattern(){
  // Constants
  double c0=3.0E8;
  // Definitions
  int Ns=1000;
  double GHz=1.0E9;
  double cm=1.0E-2;
  double theta_min=-180.0;
  double theta_max=+180.0;
  double phi_min=0.0;
  double phi_max=+360.0;
  double freq=1.05*GHz;
  double L=10*cm;
  double S=3*cm;
  double a=0.5*(1.6E-3);
  double segPerLambda=31.0;
  setQuadL(32);
  //
  Timer T=DefaultTimer;
  ticTimer(&T);
  double lambda=c0/freq;
  double theta, dtheta;
  double phi, dphi;
  dtheta = (theta_max-theta_min)/(Ns-1.0);
  dphi = (phi_max-phi_min)/(Ns-1.0);
  FILE *file=fopen("Data/BentWire/Directivity.dat", "w");
  // Choose scenario
  // createBentWire(L, S, lambda/segPerLambda); // Scenario I
  createBentWireSplit(L, S, lambda/segPerLambda); // Scenario II
  //
  Shape myShape=createShape();
  myShape.a = a;
  // Add ports
  PortsConfig myPorts=DefaultPortsConfig;
	createPortsConfig(&myPorts, 1, 0);
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
  MM_DeltaGap(&Zmn, &Vm, &In, myShape, lambda, myPorts);
  GEPP(Zmn, Vm, &In);
  double Den=computeDen(myShape, In, lambda);
  double DTheta, DPhi;
  for (int i=0; i<Ns; i++){
    phi = 0.0;
    theta = theta_min+i*dtheta;
    DTheta = DirectivityTheta(myShape, theta, phi, In, lambda, Den);
    DPhi = DirectivityPhi(myShape, theta, phi, In, lambda, Den);
    fprintf(file, "%21.14E %21.14E %21.14E %21.14E ", theta, phi, DTheta, DPhi);
    theta = 90.0;
    phi = phi_min+i*dphi;
    DTheta = DirectivityTheta(myShape, theta, phi, In, lambda, Den);
    DPhi = DirectivityPhi(myShape, theta, phi, In, lambda, Den);
    fprintf(file, "%21.14E %21.14E %21.14E %21.14E ", theta, phi, DTheta, DPhi);
    fprintf(file, "\n");
    //
    progressBar(i, Ns);
  }
  fclose(file);
  //
  deletePortsConfig(&myPorts);
  deallocateMatrix(&Zmn);
  deallocateMatrix(&Vm);
  deallocateMatrix(&In);
  deleteShape(&myShape);
  //
  unsetQuadL();
  tocTimer(&T);
  showTimer(T);
}

void programInvertedFRadiationPattern(){
  // Constants
  double c0=3.0E8;
  // Definitions
  int Ns=1000;
  double GHz=1.0E9;
  double cm=1.0E-2;
  double mm=1.0E-3;
  double theta_min=-180.0;
  double theta_max=+180.0;
  double phi_min=0.0;
  double phi_max=+360.0;
  double freq=0.925*GHz;
  double L=5.35*cm;
  double S=1.0*cm;
  double h=3.24*cm;
  double a=1.0*mm;
  double segPerLambda=51.0;
  setQuadL(32);
  //
  Timer T=DefaultTimer;
  ticTimer(&T);
  double lambda=c0/freq;
  double theta, dtheta;
  double phi, dphi;
  dtheta = (theta_max-theta_min)/(Ns-1.0);
  dphi = (phi_max-phi_min)/(Ns-1.0);
  FILE *file=fopen("Data/InvertedF/Directivity.dat", "w");
  // Choose scenario
  createInvertedF(L, S, h, lambda/segPerLambda);
  //
  Shape myShape=createShape();
  myShape.a = a;
  // Add ports
  PortsConfig myPorts=DefaultPortsConfig;
	createPortsConfig(&myPorts, 1, 0);
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
  MM_DeltaGap(&Zmn, &Vm, &In, myShape, lambda, myPorts);
  GEPP(Zmn, Vm, &In);
  double Den=computeDen(myShape, In, lambda);
  double DTheta, DPhi;
  for (int i=0; i<Ns; i++){
    phi = 0.0;
    theta = theta_min+i*dtheta;
    DTheta = DirectivityTheta(myShape, theta, phi, In, lambda, Den);
    DPhi = DirectivityPhi(myShape, theta, phi, In, lambda, Den);
    fprintf(file, "%21.14E %21.14E %21.14E %21.14E ", theta, phi, DTheta, DPhi);
    theta = 90.0;
    phi = phi_min+i*dphi;
    DTheta = DirectivityTheta(myShape, theta, phi, In, lambda, Den);
    DPhi = DirectivityPhi(myShape, theta, phi, In, lambda, Den);
    fprintf(file, "%21.14E %21.14E %21.14E %21.14E ", theta, phi, DTheta, DPhi);
    fprintf(file, "\n");
    //
    progressBar(i, Ns);
  }
  fclose(file);
  //
  deletePortsConfig(&myPorts);
  deallocateMatrix(&Zmn);
  deallocateMatrix(&Vm);
  deallocateMatrix(&In);
  deleteShape(&myShape);
  //
  unsetQuadL();
  tocTimer(&T);
  showTimer(T);
}

void programBiconicalDeltaGapRadiationPattern(){
  // Constants
  double c0=3.0E8;
  // Definitions
  int Ns=1000;
  double theta_min=-180.0;
  double theta_max=+180.0;
  double phi_min=0.0;
  double phi_max=+360.0;
  double MHz=1.0E6;
  double cm=1.0E-2;
  double mm=1.0E-3;
  double freq=78*MHz;
  int N=6;
  double l=60.35*cm;
  double gamma=87.0*mm;
  double h1=gamma/2.0;
  double h2=3.0*l/4.0;
  double h3=l/4.0;
  double S=sqrt(3.0)*l/4.0;
  double a=3.0*mm;
  double segPerLambda=31.0;
  setQuadL(32);
  //
  Timer T=DefaultTimer;
  ticTimer(&T);
  double lambda;
  double theta, dtheta;
  double phi, dphi;
  dtheta = (theta_max-theta_min)/(Ns-1.0);
  dphi = (phi_max-phi_min)/(Ns-1.0);
  FILE *file=fopen("Data/Biconical/Directivity.dat", "w");
	lambda = c0/freq;
	//
	createBiconicalAntenna(h1, h2, h3, S, N, lambda/segPerLambda);
	Shape myShape=createShape();
	myShape.a = a;
	// Add ports
	PortsConfig myPorts=DefaultPortsConfig;
	createPortsConfig(&myPorts, 1, 0);
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
	MM_DeltaGap(&Zmn, &Vm, &In, myShape, lambda, myPorts);
	GEPP(Zmn, Vm, &In);
	double Den=computeDen(myShape, In, lambda);
	double DTheta, DPhi;
  for (int i=0; i<Ns; i++){
    phi = 0.0;
    theta = theta_min+i*dtheta;
    DTheta = DirectivityTheta(myShape, theta, phi, In, lambda, Den);
    DPhi = DirectivityPhi(myShape, theta, phi, In, lambda, Den);
    fprintf(file, "%21.14E %21.14E %21.14E %21.14E ", theta, phi, DTheta, DPhi);
    theta = 90.0;
    phi = phi_min+i*dphi;
    DTheta = DirectivityTheta(myShape, theta, phi, In, lambda, Den);
    DPhi = DirectivityPhi(myShape, theta, phi, In, lambda, Den);
    fprintf(file, "%21.14E %21.14E %21.14E %21.14E ", theta, phi, DTheta, DPhi);
    fprintf(file, "\n");
    //
    progressBar(i, Ns);
  }
	//
	deletePortsConfig(&myPorts);
	deallocateMatrix(&Zmn);
	deallocateMatrix(&Vm);
	deallocateMatrix(&In);
	deleteShape(&myShape);
	//
	fclose(file);
  unsetQuadL();
  tocTimer(&T);
  showTimer(T);
}

void programBiconicalDeltaGapInputImpedance(){
  // Constants
  double c0=3.0E8;
  double MHz=1.0E6;
  double cm=1.0E-2;
  double mm=1.0E-3;
  // Definitions
  int Ns=21;
  double freq_min=25*MHz;
  double freq_max=300*MHz;
  int N=6;
  double l=60.35*cm;
  double gamma=87.0*mm;
  double h1=gamma/2.0;
  double h2=3.0*l/4.0;
  double h3=l/4.0;
  double S=sqrt(3.0)*l/4.0;
  double a=3.0*mm;
  double segPerLambda=31.0;
  setQuadL(32);
  //
  Timer T=DefaultTimer;
  ticTimer(&T);
  double lambda, freq, dfreq;
  dfreq = (freq_max-freq_min)/(Ns-1.0);
  FILE *file=fopen("Data/Biconical/Zin.dat", "w");
  for (int i=0; i<Ns; i++){
    freq = freq_min+i*dfreq;
    lambda = c0/freq;
  	//
  	createBiconicalAntenna(h1, h2, h3, S, N, lambda/segPerLambda);
  	Shape myShape=createShape();
  	myShape.a = a;
  	// Add ports
  	PortsConfig myPorts=DefaultPortsConfig;
  	createPortsConfig(&myPorts, 1, 0);
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
  	MM_DeltaGap(&Zmn, &Vm, &In, myShape, lambda, myPorts);
  	GEPP(Zmn, Vm, &In);
  	complex double Zin=0.0;
  	int flag=0;
  	for (int m=0; m<NBasis; m++){
  	  if (myShape.BasisList[m].port==1){
  		Zin=1.0/In.data[m][0];
  		flag++;
  	  }
  	}
  	assert(flag==1);
  	fprintf(file, "%21.14E %21.14E %21.14E\n", freq, creal(Zin), cimag(Zin));
    progressBar(i, Ns);
  	//
  	deletePortsConfig(&myPorts);
  	deallocateMatrix(&Zmn);
  	deallocateMatrix(&Vm);
  	deallocateMatrix(&In);
  	deleteShape(&myShape);
  }
	//
	fclose(file);
  unsetQuadL();
  tocTimer(&T);
  showTimer(T);
}
