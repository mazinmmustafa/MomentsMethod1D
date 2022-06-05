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
	
	int Ns=1000;
	double L_min, L_max, a;
	L_min = 1.0E-2;
	L_max = 4.0;
	int dL=15;
	double theta_i, phi_i;
	theta_i = 30.0;
	phi_i = 0.0;
	double E_TM, E_TE;
	E_TM = 1.0;
	E_TE = 0.0;
	
	// Correct angles
	theta_i = deg2rad(theta_i);
	phi_i = deg2rad(phi_i);
	
	// RCS
	double theta=30.0, phi=0.0;
	double sigma_theta, sigma_phi;
	
	theta = deg2rad(theta);
	phi = deg2rad(phi);
	
	double d_L=(L_max-L_min)/(Ns-1.0), L;
	
	// Save RCS data
	FILE* fID=NULL;
	char fileName[]="DataRCS.dat";
	fID = fopen(fileName, "w");
	int count=0;
	for (int i=0; i<Ns; i++){
		L = L_min+i*d_L;
		int N=round(L*dL);
		a=L/200.0;
		segment* shape=VerticalDipole(N, L); 
		TYPE** Zmn=createMatrix(N);
		TYPE* Vm=createVector(N);
		TYPE* In=createVector(N);
		MM_PlaneWave(N, a, E_TM, E_TE, theta_i, phi_i, Zmn, Vm, shape);
		GaussPivot(N, Zmn, Vm, In);
	
		sigma(N, shape, theta, phi, In, &sigma_theta, &sigma_phi);
		if (sigma_theta==0.0){
			sigma_theta = DBL_EPSILON;
		}
		if (sigma_phi==0.0){
			sigma_phi = DBL_EPSILON;
		}
		fprintf(fID, "%12.4E\t%12.4E\t%12.4E\n", L, 
		               10*log10(sigma_theta), 10*log10(sigma_phi));
		free(shape);
		deleteMatrix(N, Zmn);
		deleteVector(Vm);
		deleteVector(In);
		system("cls");
		count++;
		printf("Progress: %3d%%\n", 100*count/Ns);
	}
	fclose(fID);
	
	toc();
	
	system("Plot.gp");
}


