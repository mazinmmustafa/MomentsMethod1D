#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
// My header files
#include "Type.h"

// Matrix definitions

void zerosMatrix(int N, TYPE** A){
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			A[i][j] = (TYPE) 0;
		}
	}
};

TYPE** createMatrix(int N){
	TYPE** A=NULL;
	A = malloc(N*sizeof(TYPE*));
	for (int i=0; i<N; i++){
		A[i] = malloc(N*sizeof(TYPE));
	}
	zerosMatrix(N, A);
	return A;
}

void deleteMatrix(int N, TYPE** A){
	for (int i=0; i<N; i++){
		free(A[i]);
	}
	free(A);
}

void printMatrix(int N, TYPE** A){
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			printf("(%12.4E, %12.4E) ", creal(A[i][j]), cimag(A[i][j]));
		}
		printf("\n");
	}
}

void copyMatrix(int N, TYPE** A, TYPE** B){
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			A[i][j] = B[i][j];
		}
	}
}

// Vector Definitions

void zerosVector(int N, TYPE* V){
	for (int i=0; i<N; i++){
		V[i] = (TYPE) 0;
	}
}

TYPE* createVector(int N){
	TYPE* V=NULL;
	V = malloc(N*sizeof(TYPE));
	zerosVector(N, V);
	return V;
}

void deleteVector(TYPE* V){
	free(V);
}

void printVector(int N, TYPE* V){
	for (int i=0; i<N; i++){
		printf("(%12.4E, %12.4E)\n", creal(V[i]), cimag(V[i]));
	}
}

void copyVector(int N, TYPE* Va, TYPE* Vb){
	for (int i=0; i<N; i++){
		Va[i] = Vb[i];
	}
}