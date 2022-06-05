#ifndef SHAPE_H
#define SHAPE_H

//
#include <complex.h>
#include "../include/Vector.h"

// Definitions
typedef struct Basis{
  Vector rn_m, rn, rn_p;
  Vector Ln_m, Ln_p;
  double ln_m, ln_p;
  int port;
}Basis;
#define DefaultBasis {DefaultVector, DefaultVector, \
      DefaultVector, DefaultVector, DefaultVector, 0.0, 0.0, 0}

typedef struct Shape{
  int NBasis;
  Basis *BasisList;
  double a;
}Shape;
#define DefaultShape {0, NULL, 0.0}

typedef struct PortsConfig{
  int NExcitations, NLoads;
  int *ExcitationPorts, *LoadPorts;
  complex double *Excitations, *Loads;
}PortsConfig;
#define DefaultPortsConfig {0, 0, NULL, NULL, NULL, NULL}

// Functions
void createDipole_z(double L, double delta);
void createCircle_xy(double radius, double delta);
void createTL(double L, double S, double delta);
void createGullAntenna(double h1, double h2, double h3,
               double alpha, double delta);
void createUdaYagiAntenna(double L1, double L2, double L3,
              double S1, double S2, double S3, int N, double delta);
void createBentWire(double L, double S, double delta);
void createBentWireSplit(double L, double S, double delta);
void createInvertedF(double L, double S, double h, double delta);
void createBiconicalAntenna(double h1, double h2, double h3, double S,
            int N, double delta);
Shape createShape();
void deleteShape();
void logShape(const Shape myShape, char *fileName);
void createPortsConfig(PortsConfig *myPorts, int NExcitations, int NLoads);
void deletePortsConfig(PortsConfig *myPorts);

#endif
