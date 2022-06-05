#include <stdio.h>
// My header files
#include "Type.h"
#include "Matrix.h"

int main(){
	
	int N=4;
	
	TYPE** A=createMatrix(N);
	TYPE* V=createVector(N);
	
	printMatrix(N, A);
	printVector(N, V);
	
	deleteMatrix(N, A);
	deleteVector(V);
	
}