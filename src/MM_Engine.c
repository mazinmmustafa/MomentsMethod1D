#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
//
#include "../include/MM_Engine.h"
#include "../include/TestBench.h"
#include "../include/Vector.h"
#include "../include/Utilities.h"
#include "../include/Shape.h"
#include "../include/Matrix.h"
#include "../include/QuadL.h"

// Global variable
complex double j=I;
double pi=M_PI;
double k=2.0*M_PI;
double eta=120.0*M_PI;

complex double Rmn(double alpha, double alpha_, Basis Bm, Basis Bn,
                    double a, char s1, char s2){
    int flag=0;
    Vector r1, r2, L1, L2, B1, B2;
    if (s1=='+'&&s2=='+'){
      r1 = copyVector(Bm.rn_p);
      r2 = copyVector(Bn.rn_p);
      L1 = copyVector(Bm.Ln_p);
      L2 = copyVector(Bn.Ln_p);
      B1 = scaleVector(L1, -alpha);
      B2 = scaleVector(L2, +alpha_);
      flag++;
    }else
    if (s1=='+'&&s2=='-'){
      r1 = copyVector(Bm.rn_p);
      r2 = copyVector(Bn.rn_m);
      L1 = copyVector(Bm.Ln_p);
      L2 = copyVector(Bn.Ln_m);
      B1 = scaleVector(L1, -alpha);
      B2 = scaleVector(L2, -alpha_);
      flag++;
    }else
    if (s1=='-'&&s2=='+'){
      r1 = copyVector(Bm.rn_m);
      r2 = copyVector(Bn.rn_p);
      L1 = copyVector(Bm.Ln_m);
      L2 = copyVector(Bn.Ln_p);
      B1 = scaleVector(L1, +alpha);
      B2 = scaleVector(L2, +alpha_);
      flag++;
    }else
    if (s1=='-'&&s2=='-'){
      r1 = copyVector(Bm.rn_m);
      r2 = copyVector(Bn.rn_m);
      L1 = copyVector(Bm.Ln_m);
      L2 = copyVector(Bn.Ln_m);
      B1 = scaleVector(L1, +alpha);
      B2 = scaleVector(L2, -alpha_);
      flag++;
    }
    assert(flag==1);
    Vector A=subVector(r1, r2);
    Vector B=addVector(B1, B2);
    Vector C=addVector(A, B);
    double magC=magVector(C);
    return sqrt(magC*magC+a*a);
}

complex double gmn(double alpha, double alpha_, Basis Bm, Basis Bn,
                    double a, char s1, char s2){
    complex double R=Rmn(alpha, alpha_, Bm, Bn, a, s1, s2);
    return cexp(-j*k*R)/(4.0*pi*R);
}

typedef struct IArgs{
  Basis Bm;
  Basis Bn;
  double a;
  char s1;
  char s2;
}IArgs;

complex double Psi_mn_Integrand(double alpha, double alpha_, void *ArgsIn){
  IArgs *Args=ArgsIn;
  Basis Bm=Args->Bm;
  Basis Bn=Args->Bn;
  double a=Args->a;
  char s1=Args->s1;
  char s2=Args->s2;
  Vector Lm=DefaultVector, Ln=DefaultVector;
  complex double g=gmn(alpha, alpha_, Bm, Bn, a, s1, s2);
  int flag=0;
  if (s1=='+'&&s2=='+'){
    Lm = copyVector(Bm.Ln_p);
    Ln = copyVector(Bn.Ln_p);
    flag++;
  }else
  if (s1=='+'&&s2=='-'){
    Lm = copyVector(Bm.Ln_p);
    Ln = copyVector(Bn.Ln_m);
    flag++;
  }else
  if (s1=='-'&&s2=='+'){
    Lm = copyVector(Bm.Ln_m);
    Ln = copyVector(Bn.Ln_p);
    flag++;
  }else
  if (s1=='-'&&s2=='-'){
    Lm = copyVector(Bm.Ln_m);
    Ln = copyVector(Bn.Ln_m);
    flag++;
  }
  assert(flag==1);
  return dotVector(Lm, Ln)*alpha*alpha_*g;
}

complex double Phi_mn_Integrand(double alpha, double alpha_, void *ArgsIn){
  IArgs *Args=ArgsIn;
  Basis Bm=Args->Bm;
  Basis Bn=Args->Bn;
  double a=Args->a;
  char s1=Args->s1;
  char s2=Args->s2;
  complex double g=gmn(alpha, alpha_, Bm, Bn, a, s1, s2);
  return g;
}

complex double Psi_mn(Basis Bm, Basis Bn, double a, char s1, char s2){
  IArgs Args={Bm, Bn, a, s1, s2};
  return QuadL_2D(Psi_mn_Integrand, &Args, 0.0, 1.0, 0.0, 1.0);
}

complex double Phi_mn(Basis Bm, Basis Bn, double a, char s1, char s2){
  IArgs Args={Bm, Bn, a, s1, s2};
  return QuadL_2D(Phi_mn_Integrand, &Args, 0.0, 1.0, 0.0, 1.0);
}

double sinc(double x){
  if (fabs(x)==0.0){
    return 1.0;
  }else{
    return sin(x)/x;
  }
}

typedef struct ISingArgs{
  double L, a;
}ISingArgs;

complex double Integrand_1(double alpha, double alpha_, void *ISingArgsIn){
  ISingArgs *Args=ISingArgsIn;
  double L=Args->L;
  double a=Args->a;
  double R=sqrt(L*L*(alpha-alpha_)*(alpha-alpha_)+a*a);
  return -j*k*(L*L/(4.0*pi))*alpha*alpha_*cexp(-j*k*R/2.0)*sinc(k*R/2.0);
}

complex double Integrand_2(double alpha, double alpha_, void *ISingArgsIn){
  ISingArgs *Args=ISingArgsIn;
  double L=Args->L;
  double a=Args->a;
  double R=sqrt(L*L*(alpha-alpha_)*(alpha-alpha_)+a*a);
  return -j*k*(1.0/(4.0*pi))*cexp(-j*k*R/2.0)*sinc(k*R/2.0);
}

complex double Integrand_3(double alpha, double alpha_, void *ISingArgsIn){
  ISingArgs *Args=ISingArgsIn;
  double L=Args->L;
  double a=Args->a;
  double R=sqrt(L*L*(1.0-(alpha+alpha_))*(1.0-(alpha+alpha_))+a*a);
  return -j*k*(L*L/(4.0*pi))*alpha*alpha_*cexp(-j*k*R/2.0)*sinc(k*R/2.0);
}

complex double Integrand_4(double alpha, void *ISingArgsIn){
  ISingArgs *Args=ISingArgsIn;
  double L=Args->L;
  double a=Args->a;
  complex double A, B, C, D;
  A = sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a)+(1.0-alpha)*L;
  B = sqrt(alpha*alpha*L*L+a*a)-alpha*L;
  C = sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a);
  D = sqrt(alpha*alpha*L*L+a*a);
  return (L/(4.0*pi))*alpha*alpha*log(A/B)
        +(1.0/(4.0*pi))*alpha*(C-D);
}

complex double Integrand_5(double alpha, void *ISingArgsIn){
  ISingArgs *Args=ISingArgsIn;
  double L=Args->L;
  double a=Args->a;
  complex double A, B;
  A = sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a)+(1.0-alpha)*L;
  B = sqrt(alpha*alpha*L*L+a*a)-alpha*L;
  return (1.0/(4.0*pi*L))*log(A/B);
}

complex double Integrand_6(double alpha, void *ISingArgsIn){
  ISingArgs *Args=ISingArgsIn;
  double L=Args->L;
  double a=Args->a;
  complex double A, B, C, D;
  A = sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a)+(1.0-alpha)*L;
  B = sqrt(alpha*alpha*L*L+a*a)-alpha*L;
  C = sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a);
  D = sqrt(alpha*alpha*L*L+a*a);
  return (L/(4.0*pi))*alpha*(1.0-alpha)*log(A/B)
        -(1.0/(4.0*pi))*alpha*(C-D);
}

complex double I1(double L, double a){
  ISingArgs Args={L, a};
  return QuadL_2D(Integrand_1, &Args, 0.0, 1.0, 0.0, 1.0)
        +QuadL_1D(Integrand_4, &Args, 0.0, 1.0);
}

complex double I2(double L, double a){
  ISingArgs Args={L, a};
  return QuadL_2D(Integrand_2, &Args, 0.0, 1.0, 0.0, 1.0)
        +QuadL_1D(Integrand_5, &Args, 0.0, 1.0);
}

complex double I3(double L, double a){
  ISingArgs Args={L, a};
  return QuadL_2D(Integrand_3, &Args, 0.0, 1.0, 0.0, 1.0)
        +QuadL_1D(Integrand_6, &Args, 0.0, 1.0);
}

void MM_DeltaGap(Matrix *Zmn, Matrix *Vm, Matrix *In,
       Shape myShape, double lambda, PortsConfig myPorts){
  #define Zmn (*Zmn)
  #define Vm (*Vm)
  #define In (*In)
  // Assertions
  assert(myShape.BasisList!=NULL);
  assert(myShape.NBasis>0);
  assert(myShape.a>0.0);
  assert(Zmn.data!=NULL);
  assert(Vm.data!=NULL);
  assert(In.data!=NULL);
  assert(lambda>0.0);
  //
  double a=myShape.a/lambda;
  double tol=1.0E-10;
  int N=myShape.NBasis;
  complex double A, B;
  int flag;
  int countExcitation=0;
  int countLoad=0;
  for (int m=0; m<N; m++){
    // Obtain Bm
    Basis Bm={scaleVector(myShape.BasisList[m].rn_m, 1.0/lambda),
              scaleVector(myShape.BasisList[m].rn, 1.0/lambda),
              scaleVector(myShape.BasisList[m].rn_p, 1.0/lambda),
              scaleVector(myShape.BasisList[m].Ln_m, 1.0/lambda),
              scaleVector(myShape.BasisList[m].Ln_p, 1.0/lambda),
              myShape.BasisList[m].ln_m/lambda,
              myShape.BasisList[m].ln_p/lambda,
              myShape.BasisList[m].port};
    // Add ports
    if (countExcitation<myPorts.NExcitations){
      if (Bm.port==myPorts.ExcitationPorts[countExcitation]){
        Vm.data[m][0] = myPorts.Excitations[countExcitation];
        countExcitation++;
      }
    }
    for (int n=0; n<N; n++){
      // Obtain Bn
      Basis Bn={scaleVector(myShape.BasisList[n].rn_m, 1.0/lambda),
                scaleVector(myShape.BasisList[n].rn, 1.0/lambda),
                scaleVector(myShape.BasisList[n].rn_p, 1.0/lambda),
                scaleVector(myShape.BasisList[n].Ln_m, 1.0/lambda),
                scaleVector(myShape.BasisList[n].Ln_p, 1.0/lambda),
                myShape.BasisList[n].ln_m/lambda,
                myShape.BasisList[n].ln_p/lambda,
                myShape.BasisList[n].port};
      // Check scenarios
      flag = 0;
      if (isVectorEqual(Bm.rn_m, Bn.rn_m, tol)&&
          isVectorEqual(Bm.rn, Bn.rn, tol)&&
          isVectorEqual(Bm.rn_p, Bn.rn_p, tol)){
        A = I1(Bm.ln_p, a)
           +Psi_mn(Bm, Bn, a, '+', '-')
           +Psi_mn(Bm, Bn, a, '-', '+')
           +I1(Bm.ln_m, a);

        B = I2(Bm.ln_p, a)
           -Phi_mn(Bm, Bn, a, '+', '-')
           -Phi_mn(Bm, Bn, a, '-', '+')
           +I2(Bm.ln_m, a);
        flag++;
        // printf("Scenario I\n");
      }else
      if (isVectorEqual(Bm.rn, Bn.rn_m, tol)&&
          isVectorEqual(Bm.rn_p, Bn.rn, tol)){
        A = Psi_mn(Bm, Bn, a, '+', '+')
           +I3(Bm.ln_p, a)
           +Psi_mn(Bm, Bn, a, '-', '+')
           +Psi_mn(Bm, Bn, a, '-', '-');

        B = Phi_mn(Bm, Bn, a, '+', '+')
           -I2(Bm.ln_p, a)
           -Phi_mn(Bm, Bn, a, '-', '+')
           +Phi_mn(Bm, Bn, a, '-', '-');
        flag++;
        // printf("Scenario II\n");
      }else
      if (isVectorEqual(Bm.rn_m, Bn.rn, tol)&&
          isVectorEqual(Bm.rn, Bn.rn_p, tol)){
        A = Psi_mn(Bm, Bn, a, '+', '+')
           +Psi_mn(Bm, Bn, a, '+', '-')
           +I3(Bm.ln_m, a)
           +Psi_mn(Bm, Bn, a, '-', '-');

        B = Phi_mn(Bm, Bn, a, '+', '+')
           -Phi_mn(Bm, Bn, a, '+', '-')
           -I2(Bm.ln_m, a)
           +Phi_mn(Bm, Bn, a, '-', '-');
        flag++;
        // printf("Scenario III\n");
      }else{
        A = Psi_mn(Bm, Bn, a, '+', '+')
           +Psi_mn(Bm, Bn, a, '+', '-')
           +Psi_mn(Bm, Bn, a, '-', '+')
           +Psi_mn(Bm, Bn, a, '-', '-');

        B = Phi_mn(Bm, Bn, a, '+', '+')
           -Phi_mn(Bm, Bn, a, '+', '-')
           -Phi_mn(Bm, Bn, a, '-', '+')
           +Phi_mn(Bm, Bn, a, '-', '-');
        flag++;
        // printf("Scenario IV\n");
      }
      assert(flag==1);
      Zmn.data[m][n] = j*k*eta*A-j*(eta/k)*B;
      // Add loads
      if (countLoad<myPorts.NLoads&&m==n){
        if (Bm.port==myPorts.LoadPorts[countLoad]){
          Zmn.data[m][n]+=myPorts.Loads[countLoad];
          countLoad++;
        }
      }
    }
  }
  //
  #undef Zmn
  #undef Vm
  #undef In
}


