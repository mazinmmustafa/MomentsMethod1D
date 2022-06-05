#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
// My header files
#include "Type.h"
#include "Matrix.h"
#include "GEEP.h"

int main(){
	
	int N=3;
	
	TYPE** A=createMatrix(N);
	TYPE* x=createVector(N);
	TYPE* b=createVector(N);
	
	A[0][0] = 1.0+I*4.0; A[0][1] = 1.0+I*2.0; A[0][2] = 1.0+I*2.0;
	A[1][0] = 1.0+I*4.0; A[1][1] = 1.0+I*3.0; A[1][2] = 1.0+I*4.0;
	A[2][0] = 1.0+I*3.0; A[2][1] = 1.0+I*4.0; A[2][2] = 1.0+I*3.0;
	
	printf("A = \n");
	printMatrix(N, A);
	
	b[0] = 1.0+I*2.0; 
	b[1] = 2.0+I*4.0; 
	b[2] = 1.0+I*4.0;
	
	printf("b = \n");
	printVector(N, b);
	
	GaussPivot(N, A, b, x);
	
	printf("x = \n");
	printVector(N, x);
	
	deleteMatrix(N, A);
	deleteVector(x);
	deleteVector(b);
	
}