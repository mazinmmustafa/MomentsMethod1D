#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

TYPE ZinVerticalDipole(int N, double L, double a){
	int M;
	segment* shape=VerticalDipoleDeltaGap(&N, L, &M); 
	TYPE** Zmn=createMatrix(N);
	TYPE* Vm=createVector(N);
	TYPE* In=createVector(N);
	MM_DeltaGap(N, M, a, Zmn, Vm, shape);
	GaussPivot(N, Zmn, Vm, In);
	TYPE Zin=1.0/In[M];
	free(shape);
	deleteMatrix(N, Zmn);
	deleteVector(Vm);
	deleteVector(In);
	return Zin;
};

int main(){
	
	tic();
	
	int Ns;
	double L_min, L_max;
	double ratio;

	Ns = 601;
	L_min = 1.0E-2;
	L_max = 2.0E+0;
	
	ratio = 2.0*74.2;
	
	int N_per_lambda=15;
	
	FILE* fID=NULL;
	char fileName[]="Data.dat";
	fID = fopen(fileName, "w");
	
	TYPE Yin;
	int N;
	double a, L, dL=(L_max-L_min)/(Ns-1.0);
	
	for (int i=0; i<Ns; i++){
		L = L_min+i*dL;
		a = L/ratio;
		N = round(N_per_lambda*L);
		if (N==0){
			N++;
		}
		Yin = 1E3/ZinVerticalDipole(N, L, a);
		fprintf(fID, "%12.4E\t%12.4E\t%12.4E\n", L, creal(Yin), cimag(Yin));
	}
	
	fclose(fID);
	
	toc();
	
}


