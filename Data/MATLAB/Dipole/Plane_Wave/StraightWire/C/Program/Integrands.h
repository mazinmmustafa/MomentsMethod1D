#ifndef INTEGRANDS_H
#define INTEGRANDS_H

#include "Type.h"
#include "Segment.h"

TYPE psi_pp(segment, segment, double, double);
TYPE psi_pm(segment, segment, double, double);
TYPE psi_mp(segment, segment, double, double);
TYPE psi_mm(segment, segment, double, double);

TYPE phi_pp(segment, segment, double, double);
TYPE phi_pm(segment, segment, double, double);
TYPE phi_mp(segment, segment, double, double);
TYPE phi_mm(segment, segment, double, double);

#endif