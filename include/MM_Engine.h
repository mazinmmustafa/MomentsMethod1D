#ifndef MM_ENGINE_H
#define MM_ENGINE_H

//
#include "../include/Matrix.h"
#include "../include/Shape.h"

// Functions
void MM_DeltaGap(Matrix *Zmn, Matrix *Vm, Matrix *In,
         Shape myShape, double lambda, PortsConfig myPorts);

#endif