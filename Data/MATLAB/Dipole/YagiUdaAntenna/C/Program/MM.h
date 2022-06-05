#ifndef MM_H
#define MM_H

#include "Type.h"

void MM_PlaneWave(int, double, TYPE, TYPE, double, double, TYPE**, TYPE*, const segment*);
void MM_DeltaGap(int, int, double, TYPE**, TYPE*, const segment*);

#endif