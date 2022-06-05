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
#include "AMath.h"
#include "RCS.h"

int main(){
	
	tic();
	
	double L=3.0;
	double L1=0.7*L, L2=0.3*L, a=L/200.0;
	int dL=15;
	double theta_i, phi_i;
	theta_i = 30.0;
	phi_i = 30.0;
	double E_TM, E_TE;
	E_TM = 1.0;
	E_TE = 0.0;
	
	int N1=round(L1*dL), N2=round(L2*dL);
	
	// Correct angles
	theta_i = deg2rad(theta_i);
	phi_i = deg2rad(phi_i);
	
	segment* shape=LShapeWire(N1, N2, L1, L2); 
	int N=N1+N2+1;
	
	// printShape(N, shape);
	
	TYPE** Zmn=createMatrix(N);
	TYPE* Vm=createVector(N);
	TYPE* In=createVector(N);

	MM_PlaneWave(N, a, E_TM, E_TE, theta_i, phi_i, Zmn, Vm, shape);
	GaussPivot(N, Zmn, Vm, In);
	
	// Save currents data
	FILE* fID1=NULL;
	char fileName1[]="DataCurrent.dat";
	fID1 = fopen(fileName1, "w");
	for (int n=0; n<N; n++){
		fprintf(fID1, "%12.4E\n", cabs(In[n]));
	}
	fclose(fID1);

	// RCS
	int Ns=1000;
	double phi_min, phi_max, d_phi, phi;
	phi_min = -90.0;
	phi_max = +90.0;
	double theta=60.0;
	
	theta = deg2rad(theta);
	phi_min = deg2rad(phi_min);
	phi_max = deg2rad(phi_max);
	d_phi = (phi_max-phi_min)/(Ns-1.0);
	double sigma_theta, sigma_phi;
	
	// Save RCS data
	FILE* fID2=NULL;
	char fileName2[]="DataRCS.dat";
	fID2 = fopen(fileName2, "w");
	for (int i=0; i<Ns; i++){
		phi = phi_min+i*d_phi;
		sigma(N, shape, theta, phi, In, &sigma_theta, &sigma_phi);
		if (sigma_theta==0.0){
			sigma_theta = DBL_EPSILON;
		}
		if (sigma_phi==0.0){
			sigma_phi = DBL_EPSILON;
		}
		fprintf(fID2, "%12.4E\t%12.4E\t%12.4E\n", rad2deg(phi), 
		               10*log10(sigma_theta), 10*log10(sigma_phi));
	}
	fclose(fID2);
	
	free(shape);
	deleteMatrix(N, Zmn);
	deleteVector(Vm);
	deleteVector(In);
	
	toc();

	system("Plot.gp");
	
}


