#include <stdio.h>
#include <math.h>
#include <complex.h>
//
#include "../include/Vector.h"

void showVector(Vector A){
  printf("[(%9.2E, %9.2E), (%9.2E, %9.2E), (%9.2E, %9.2E)]\n",
        creal(A.x), cimag(A.x),
        creal(A.y), cimag(A.y),
        creal(A.z), cimag(A.z));
}

void showRealVector(RealVector A){
  printf("[%9.2E, %9.2E, %9.2E]\n",
          A.x, A.y, A.z);
}

double magVector(Vector A){
  return sqrt(cabs(A.x*A.x)+
              cabs(A.y*A.y)+
              cabs(A.z*A.z));
}

double magRealVector(RealVector A){
  return sqrt(A.x*A.x+A.y*A.y+A.z*A.z);
}

complex double dotVector(Vector A, Vector B){
  return A.x*B.x+A.y*B.y+A.z*B.z;
}

double dotRealVector(RealVector A, RealVector B){
  return A.x*B.x+A.y*B.y+A.z*B.z;
}

Vector crossVector(Vector A, Vector B){
  Vector C={A.y*B.z-A.z*B.y,
            A.z*B.x-A.x*B.z,
            A.x*B.y-A.y*B.x};
  return C;
}

RealVector crossRealVector(RealVector A, RealVector B){
  RealVector C={A.y*B.z-A.z*B.y,
                A.z*B.x-A.x*B.z,
                A.x*B.y-A.y*B.x};
  return C;
}

RealVector unitRealVector(RealVector A){
  double magA=magRealVector(A);
  RealVector n={A.x/magA, A.y/magA, A.z/magA};
  return n;
}

Vector addVector(Vector A, Vector B){
  Vector C={A.x+B.x, A.y+B.y, A.z+B.z};
  return C;
}

RealVector addRealVector(RealVector A, RealVector B){
  RealVector C={A.x+B.x, A.y+B.y, A.z+B.z};
  return C;
}

Vector subVector(Vector A, Vector B){
  Vector C={A.x-B.x, A.y-B.y, A.z-B.z};
  return C;
}

RealVector subRealVector(RealVector A, RealVector B){
  RealVector C={A.x-B.x, A.y-B.y, A.z-B.z};
  return C;
}

Vector scaleVector(Vector A, complex double a){
  Vector B={a*A.x, a*A.y, a*A.z};
  return B;
}

RealVector scaleRealVector(RealVector A, double a){
  RealVector B={a*A.x, a*A.y, a*A.z};
  return B;
}

Vector RealVector2Vector(RealVector A){
  Vector B={A.x, A.y, A.z};
  return B;
}

RealVector Vector2RealVector(Vector A){
  RealVector B={creal(A.x), creal(A.y), creal(A.z)};
  return B;
}

int isVectorEqual(Vector A, Vector B, double tol){
  if (magVector(subVector(A, B))<tol){
    return 1;
  }else{
    return 0;
  }
}

int isRealVectorEqual(RealVector A, RealVector B, double tol){
  if (magRealVector(subRealVector(A, B))<tol){
    return 1;
  }else{
    return 0;
  }
}

Vector copyVector(Vector A){
  Vector B={A.x, A.y, A.z};
  return B;
}

RealVector copyRealVector(RealVector A){
  RealVector B={A.x, A.y, A.z};
  return B;
}

