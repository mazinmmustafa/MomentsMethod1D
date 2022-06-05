#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
//
#include "../include/FarField.h"
#include "../include/Vector.h"
#include "../include/Utilities.h"
#include "../include/Shape.h"
#include "../include/Matrix.h"
#include "../include/GaussLegendreRule.h"

complex double IRadiation(double a, double b, char s){
  complex double j=I;
  complex double ans=0.0;
  if (fabs(a)<1.0E-4){
    ans = cexp(j*b)/2.0;
  }else{
    if (s=='+'){
      ans = cexp(+j*b)*(cexp(-j*a)*(1.0+j*a)-1.0)/(a*a);
    }
    if (s=='-'){
      ans = cexp(+j*b)*(cexp(+j*a)*(1.0-j*a)-1.0)/(a*a);
    }
  }
  return ans;
}

void sigma(Shape myShape, double theta, double phi,
          const Matrix In, double lambda,
          double *sigmaTheta, double *sigmaPhi){
  // Assertions
  assert(myShape.NBasis>0);
  assert(myShape.BasisList!=NULL);
  // Constants
  double pi=M_PI;
  double k=2.0*pi;
  double eta=120.0*pi;
  //
  theta = deg2rad(theta);
  phi = deg2rad(phi);
  //
  Vector r_hat={sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
  Vector theta_hat={cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)};
  Vector phi_hat={-sin(phi), cos(phi), 0.0};
  //
  complex double sumTheta=0.0, sumPhi=0.0;
  complex double kn_m, kn_p;
  for (int n=0; n<myShape.NBasis; n++){
    Basis Bn={scaleVector(myShape.BasisList[n].rn_m, 1.0/lambda),
              scaleVector(myShape.BasisList[n].rn, 1.0/lambda),
              scaleVector(myShape.BasisList[n].rn_p, 1.0/lambda),
              scaleVector(myShape.BasisList[n].Ln_m, 1.0/lambda),
              scaleVector(myShape.BasisList[n].Ln_p, 1.0/lambda),
              myShape.BasisList[n].ln_m/lambda,
              myShape.BasisList[n].ln_p/lambda,
              myShape.BasisList[n].port};
    kn_p = IRadiation(k*dotVector(Bn.Ln_p, r_hat), k*dotVector(Bn.rn_p, r_hat), '+');
    kn_m = IRadiation(k*dotVector(Bn.Ln_m, r_hat), k*dotVector(Bn.rn_m, r_hat), '-');
    sumTheta+=In.data[n][0]*(dotVector(theta_hat, scaleVector(Bn.Ln_p, kn_p))
                            +dotVector(theta_hat, scaleVector(Bn.Ln_m, kn_m)));
    sumPhi+=In.data[n][0]*(dotVector(phi_hat, scaleVector(Bn.Ln_p, kn_p))
                          +dotVector(phi_hat, scaleVector(Bn.Ln_m, kn_m)));
  }
  (*sigmaTheta) = pi*eta*eta*cabs(sumTheta)*cabs(sumTheta);
  (*sigmaPhi) = pi*eta*eta*cabs(sumPhi)*cabs(sumPhi);
}

typedef struct FArgs FArgs;
struct FArgs{
  Shape myShape;
  Matrix In;
  double lambda;
};

double FTheta(double theta, double phi, void *Args){
  //
  FArgs *myArgs=Args;
  Shape myShape=myArgs->myShape;
  Matrix In=myArgs->In;
  double lambda=myArgs->lambda;
  // Assertions
  assert(myShape.NBasis>0);
  assert(myShape.BasisList!=NULL);
  // Constants
  double pi=M_PI;
  double k=2.0*pi;
  //
  Vector r_hat={sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
  Vector theta_hat={cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)};
  //
  complex double sumTheta=0.0;
  complex double kn_m, kn_p;
  for (int n=0; n<myShape.NBasis; n++){
    Basis Bn={scaleVector(myShape.BasisList[n].rn_m, 1.0/lambda),
              scaleVector(myShape.BasisList[n].rn, 1.0/lambda),
              scaleVector(myShape.BasisList[n].rn_p, 1.0/lambda),
              scaleVector(myShape.BasisList[n].Ln_m, 1.0/lambda),
              scaleVector(myShape.BasisList[n].Ln_p, 1.0/lambda),
              myShape.BasisList[n].ln_m/lambda,
              myShape.BasisList[n].ln_p/lambda,
              myShape.BasisList[n].port};
    kn_p = IRadiation(k*dotVector(Bn.Ln_p, r_hat), k*dotVector(Bn.rn_p, r_hat), '+');
    kn_m = IRadiation(k*dotVector(Bn.Ln_m, r_hat), k*dotVector(Bn.rn_m, r_hat), '-');
    sumTheta+=In.data[n][0]*(dotVector(theta_hat, scaleVector(Bn.Ln_p, kn_p))
                            +dotVector(theta_hat, scaleVector(Bn.Ln_m, kn_m)));
  }
  return cabs(sumTheta)*cabs(sumTheta);
}

double FPhi(double theta, double phi, void *Args){
  //
  FArgs *myArgs=Args;
  Shape myShape=myArgs->myShape;
  Matrix In=myArgs->In;
  double lambda=myArgs->lambda;
  // Assertions
  assert(myShape.NBasis>0);
  assert(myShape.BasisList!=NULL);
  // Constants
  double pi=M_PI;
  double k=2.0*pi;
  //
  Vector r_hat={sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
  Vector phi_hat={-sin(phi), cos(phi), 0.0};
  //
  complex double sumPhi=0.0;
  complex double kn_m, kn_p;
  for (int n=0; n<myShape.NBasis; n++){
    Basis Bn={scaleVector(myShape.BasisList[n].rn_m, 1.0/lambda),
              scaleVector(myShape.BasisList[n].rn, 1.0/lambda),
              scaleVector(myShape.BasisList[n].rn_p, 1.0/lambda),
              scaleVector(myShape.BasisList[n].Ln_m, 1.0/lambda),
              scaleVector(myShape.BasisList[n].Ln_p, 1.0/lambda),
              myShape.BasisList[n].ln_m/lambda,
              myShape.BasisList[n].ln_p/lambda,
              myShape.BasisList[n].port};
    kn_p = IRadiation(k*dotVector(Bn.Ln_p, r_hat), k*dotVector(Bn.rn_p, r_hat), '+');
    kn_m = IRadiation(k*dotVector(Bn.Ln_m, r_hat), k*dotVector(Bn.rn_m, r_hat), '-');
    sumPhi+=In.data[n][0]*(dotVector(phi_hat, scaleVector(Bn.Ln_p, kn_p))
                          +dotVector(phi_hat, scaleVector(Bn.Ln_m, kn_m)));
  }
  return cabs(sumPhi)*cabs(sumPhi);
}

double QuadL_2D_Custom(double func(double, double, void*),
		void *funcArgs, int NQuad){
	// \int_{a2}^{b2} \int_{a1}^{b1} func(x, y) dx dy
  double *xL=malloc(NQuad*sizeof(double));
  double *wL=malloc(NQuad*sizeof(double));
  GaussLegendreRule(NQuad, xL, wL);
  // Constants
  double pi=M_PI;
  //
  double h1=pi/2.0;
	double h2=pi;
	complex double sum=0.0;
	for (int i=0; i<NQuad; i++){
		for (int j=0; j<NQuad; j++){
			sum+=h1*h2*wL[i]*wL[j]*func(h1*xL[i]+h1, h2*xL[j]+h2, funcArgs)*sin(h1*xL[i]+h1);
		}
	}
  //
  free(xL);
  free(wL);
  xL = NULL;
  wL = NULL;
	return sum;
}

double Integrand(double theta, double phi, void *Args){
  return FTheta(theta, phi, Args)+FPhi(theta, phi, Args);
}

double computeDen(Shape myShape, Matrix In, double lambda){
  //
  FArgs myArgs;
  myArgs.myShape = myShape;
  myArgs.In = In;
  myArgs.lambda = lambda;
  int NQuad=512;
  return QuadL_2D_Custom(Integrand, &myArgs, NQuad);
}

double DirectivityTheta(Shape myShape, double theta, double phi,
          const Matrix In, double lambda, double Den){
  // Constants
  double pi=M_PI;
  //
  FArgs myArgs;
  myArgs.myShape = myShape;
  myArgs.In = In;
  myArgs.lambda = lambda;
  //
  theta = deg2rad(theta);
  phi = deg2rad(phi);
  //
  return 4.0*pi*FTheta(theta, phi, &myArgs)/Den;
}

double DirectivityPhi(Shape myShape, double theta, double phi,
          const Matrix In, double lambda, double Den){
  // Constants
  double pi=M_PI;
  //
  FArgs myArgs;
  myArgs.myShape = myShape;
  myArgs.In = In;
  myArgs.lambda = lambda;
  //
  theta = deg2rad(theta);
  phi = deg2rad(phi);
  //
  return 4.0*pi*FPhi(theta, phi, &myArgs)/Den;
}
