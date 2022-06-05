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
#include "FarField.h"

int main(){
	
	tic();
	
	int Nd=13;
	double Lr=0.5, Lf=0.47, Ld=0.406, a=0.003;
	double Sr=0.25, Sd=0.34;
	double dL=31; // More than 101 is needed!
	int N, M;
	int *M_elements=malloc((Nd+2)*sizeof(int));
	
	segment* shape=YagiUdaAntenna(&N, dL, Nd, Lr, Lf, Ld, Sr, Sd, &M, M_elements); 
	// printShape(N, shape);
	printf("M = %d\n", M);
	
	TYPE** Zmn=createMatrix(N);
	TYPE* Vm=createVector(N);
	TYPE* In=createVector(N);

	MM_DeltaGap(N, M, a, Zmn, Vm, shape);
	GaussPivot(N, Zmn, Vm, In);
	
	// Save currents data
	FILE* fID1=NULL;
	char fileName1[]="DataCurrent.dat";
	fID1 = fopen(fileName1, "w");
	for (int n=0; n<N; n++){
		fprintf(fID1, "%12.4E\n", cabs(In[n])/1E-3);
	}
	fclose(fID1);
	
	double* InPeak=malloc((Nd+2)*sizeof(double));
	for (int i=0; i<(Nd+2); i++){
		InPeak[i] = cabs(In[M_elements[i]])/1E-3;
	}
	
	double Imax=InPeak[0];
	for (int i=0; i<(Nd+2); i++){
		if (InPeak[i]>Imax){
			Imax = InPeak[i];
		}
	}
	for (int i=0; i<(Nd+2); i++){
		InPeak[i] /= Imax;
	}
	
	// Save peak currents data
	FILE* fID2=NULL;
	char fileName2[]="DataPeakCurrent.dat";
	fID2 = fopen(fileName2, "w");
	for (int i=0; i<(Nd+2); i++){
		fprintf(fID2, "%12.4E\n", InPeak[i]);
	}
	fclose(fID2);

	int Ns=1000;
	double theta_min, theta_max, d_theta, theta;
	double phi_min, phi_max, d_phi, phi;
	double E_max, E_theta, E_phi;;
	// Elevation
	
	theta_min = 0.0;
	theta_max = 360.0;
	phi = 0.0;
	phi = deg2rad(phi);
	
	theta_min = deg2rad(theta_min);
	theta_max = deg2rad(theta_max);
	d_theta = (theta_max-theta_min)/(Ns-1.0);
	
	double* E_theta_Far=malloc(Ns*sizeof(double));
	for (int i=0; i<Ns; i++){
		theta = theta_min+i*d_theta;
		FarField(N, shape, theta, phi, In, &E_theta, &E_phi);
		E_theta_Far[i] = E_theta;
		if (E_theta_Far[i]==0.0){
			E_theta_Far[i] = DBL_EPSILON;
		}
		if (E_theta_Far[i]==0.0){
			E_theta_Far[i] = DBL_EPSILON;
		}
	}
	E_max=E_theta_Far[0];
	printf("%f\n", E_max);
	for (int i=0; i<Ns; i++){
		if (E_max<E_theta_Far[i]){
			E_max = E_theta_Far[i];
		}
	}
	for (int i=0; i<Ns; i++){
		E_theta_Far[i] /= E_max;
	}
	
	// Azimuth
	phi_min = 0.0;
	phi_max = 360.0;
	theta = 90.0;
	theta = deg2rad(theta);
	
	phi_min = deg2rad(phi_min);
	phi_max = deg2rad(phi_max);
	d_phi = (phi_max-phi_min)/(Ns-1.0);
	
	double* E_phi_Far=malloc(Ns*sizeof(double));
	for (int i=0; i<Ns; i++){
		phi = phi_min+i*d_phi;
		FarField(N, shape, theta, phi, In, &E_theta, &E_phi);
		E_phi_Far[i] = E_theta;
		if (E_phi_Far[i]==0.0){
			E_phi_Far[i] = DBL_EPSILON;
		}
		if (E_phi_Far[i]==0.0){
			E_phi_Far[i] = DBL_EPSILON;
		}
	}
	E_max=E_phi_Far[0];
	printf("%f\n", E_max);
	for (int i=0; i<Ns; i++){
		if (E_max<E_phi_Far[i]){
			E_max = E_phi_Far[i];
		}
	}
	for (int i=0; i<Ns; i++){
		E_phi_Far[i] /= E_max;
	}
	
	// Save far field data
	FILE* fID3=NULL;
	char fileName3[]="DataElevation.dat";
	fID3 = fopen(fileName3, "w");
	for (int i=0; i<Ns; i++){
		theta = theta_min+i*d_theta;
		phi = phi_min+i*d_phi;
		fprintf(fID3, "%12.4E\t%12.4E%12.4E\t%12.4E\n", 
		         theta, E_theta_Far[i], phi, E_phi_Far[i]);
	}
	fclose(fID3);
	
	free(shape);
	deleteMatrix(N, Zmn);
	deleteVector(Vm);
	deleteVector(In);
	free(M_elements);
	free(InPeak);
	toc();

	system("Plot.gp");
	
}


