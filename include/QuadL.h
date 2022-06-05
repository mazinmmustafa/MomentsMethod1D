#ifndef QUADL_H
#define QUADL_H

//
#include <complex.h>

// Functions
void setQuadL(int N);
void unsetQuadL();
complex double QuadL_1D(complex double func(double, void*), 
        void *funcArgs, double a, double b);
complex double QuadL_2D(complex double func(double, double, void*), 
        void *funcArgs, double a2, double b2, double a1, double b1);
complex double QuadL_Tri_2D(complex double func(double, double, void*), void *funcArgs);
complex double QuadL_Tri_4D(complex double func(double, double, double, double, void*), 
        void *funcArgs);

#endif