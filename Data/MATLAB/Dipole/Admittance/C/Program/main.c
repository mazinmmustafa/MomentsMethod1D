#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
// My header files
#include "Type.h"
#include "Segment.h"
#include "Shapes.h"
#include "Matrix.h"
#include "MM.h"
#include "GEEP.h"
#include "Timer.h"

int main(){
	
	tic();
	
	double L=0.47, a=0.005;
	int N=31, M;
	
	// User input
	printf("Enter L: ");
	scanf("%lf", &L);
	printf("Enter a: ");
	scanf("%lf", &a);
	printf("Enter N: ");
	scanf("%d", &N);
	
	// Test vertical dipole shape
	segment* shape=VerticalDipoleDeltaGap(&N, L, &M); 
	// printShape(N, shape); // Only if you want to print
	
	TYPE** Zmn=createMatrix(N);
	TYPE* Vm=createVector(N);
	TYPE* In=createVector(N);
	
	// Test Method of Moments
	MM_DeltaGap(N, M, a, Zmn, Vm, shape);
	
	// printf("\nZ = \n");
	// printMatrix(N, Zmn); // Only if you want to print
	// printf("\nV = \n");
	// printVector(N, Vm); // Only if you want to print
	
	// Test currents
	GaussPivot(N, Zmn, Vm, In);
	
	// Save currents data
	FILE* fID=NULL;
	char fileName[]="Data.dat";
	fID = fopen(fileName, "w");
	for (int n=0; n<N; n++){
		fprintf(fID, "%12.4E\t%12.4E\n", creal(In[n]), cimag(In[n]));
	}
	fclose(fID);
	
	// printf("\nI = \n");
	// printVector(N, In); // Only if you want to print
	
	TYPE Zin=1.0/In[M];
	
	printf("\nZin = (%6.2f, %6.2f) [ohm]\n", creal(Zin), cimag(Zin));
	
	free(shape);
	deleteMatrix(N, Zmn);
	deleteVector(Vm);
	deleteVector(In);
	
	toc();

}


