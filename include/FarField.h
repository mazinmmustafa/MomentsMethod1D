#ifndef FARFIELD_H
#define FARFIELD_H

//
#include <complex.h>
//
#include "../include/Shape.h"
#include "../include/Matrix.h"

// Functions
complex double IRadiation(double a, double b, char s);
void sigma(Shape myShape, double theta, double phi,
          const Matrix In, double lambda,
          double *sigmaTheta, double *sigmaPhi);
double computeDen(Shape myShape, Matrix In, double lambda);
double DirectivityTheta(Shape myShape, double theta, double phi,
          const Matrix In, double lambda, double Den);
double DirectivityPhi(Shape myShape, double theta, double phi,
          const Matrix In, double lambda, double Den);
              
#endif
