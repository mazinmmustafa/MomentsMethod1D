#ifndef MATRIX_H
#define MATRIX_H

#include "Type.h"

// Matrix
void zerosMatrix(int, TYPE**);
TYPE** createMatrix(int);
void deleteMatrix(int, TYPE**);
void printMatrix(int, TYPE**);
void copyMatrix(int, TYPE**, TYPE**);

// Vector
void zerosVector(int, TYPE*);
TYPE* createVector(int);
void deleteVector(TYPE*);
void printVector(int, TYPE*);
void copyVector(int, TYPE*, TYPE*);

#endif