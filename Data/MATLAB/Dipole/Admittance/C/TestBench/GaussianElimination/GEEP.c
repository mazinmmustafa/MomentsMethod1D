#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
// My header files
#include "Type.h"
#include "Matrix.h"

void PivotSwap(int N, TYPE** A, TYPE* b, int k){
	TYPE Amax;
	int p;
	Amax = A[k][k];
	p = k;
	for (int i=k+1; i<N; i++){
		if (cabs(A[i][k])>cabs(Amax)){
			Amax = A[i][k];
			p = i;
		}
	}
	TYPE D;
	if (p!=k){
		D = b[p];
		b[p] = b[k];
		b[k] = D;
		for (int j=0; j<N; j++){
			D = A[p][j];
			A[p][j] = A[k][j];
			A[k][j] = D;
		}
	}
}

int GaussPivot(int N, TYPE** A_in, TYPE* b_in, TYPE* x){
	TYPE** A=createMatrix(N);
	TYPE* b=createVector(N);
	copyMatrix(N, A, A_in);
	copyVector(N, b, b_in);
	for (int k=0; k<N-1; k++){
		PivotSwap(N, A, b, k);
		if (cabs(A[k][k])==0.0){
			printf("No solution found!\n");
			return -1;
		}
		for (int i=k+1; i<N; i++){
			b[i] -= A[i][k]*b[k]/A[k][k];
			for (int j=k+1; j<N; j++){
				A[i][j] -= A[i][k]*A[k][j]/A[k][k];
			}
			A[i][k] = (TYPE) 0.0;
		}
	}
	TYPE s;
	for (int i=N-1; i>-1; i--){
		s = (TYPE) 0.0;
		for (int j=i+1; j<N; j++){
			s += A[i][j]*x[j];
		}
		x[i] = (b[i]-s)/A[i][i];
	}
	return 0;
}