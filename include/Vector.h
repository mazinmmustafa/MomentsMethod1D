#ifndef VECTOR_H
#define VECTOR_H

// Includes
#include <complex.h>

// Definitions
typedef struct Vector{
  complex double x, y, z;
}Vector;
#define DefaultVector {0.0, 0.0, 0.0}

typedef struct RealVector{
  double x, y, z;
}RealVector;
#define DefaultRealVector {0.0, 0.0, 0.0}

// Functions
void showVector(Vector A);
void showRealVector(RealVector A);
double magVector(Vector A);
double magRealVector(RealVector A);
complex double dotVector(Vector A, Vector B);
double dotRealVector(RealVector A, RealVector B);
Vector crossVector(Vector A, Vector B);
RealVector crossRealVector(RealVector A, RealVector B);
RealVector unitRealVector(RealVector A);
Vector addVector(Vector A, Vector B);
RealVector addRealVector(RealVector A, RealVector B);
Vector subVector(Vector A, Vector B);
RealVector subRealVector(RealVector A, RealVector B);
Vector scaleVector(Vector A, complex double a);
RealVector scaleRealVector(RealVector A, double a);
Vector RealVector2Vector(RealVector A);
RealVector Vector2RealVector(Vector A);
int isVectorEqual(Vector A, Vector B, double tol);
int isRealVectorEqual(RealVector A, RealVector B, double tol);
Vector copyVector(Vector A);
RealVector copyRealVector(RealVector A);

#endif
