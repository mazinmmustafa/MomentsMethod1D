#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// My header files
#include "Segment.h"

segment* YagiUdaAntenna(int* N, double dL_est, int Nd, double Lr, double Lf, 
                   double Ld, double Sr, double Sd, int* M, int* M_elements){
	int N_elements[Nd+2];
	int intSum=0;
	N_elements[0]= round(Lr*dL_est);
	if ((N_elements[0])%2==0){(N_elements[0])++;}
	(M_elements[0])=((N_elements[0])-1)/2;
	N_elements[1]= round(Lf*dL_est);
	if ((N_elements[1])%2==0){(N_elements[1])++;}
	intSum += N_elements[0];
	(M_elements[1])=intSum+((N_elements[1])-1)/2;
	for (int i=2; i<(Nd+2); i++){
		N_elements[i] = round(Ld*dL_est);
		if ((N_elements[i])%2==0){(N_elements[i])++;}
		intSum += N_elements[i-1];
		(M_elements[i])=intSum+((N_elements[i])-1)/2;
	}
	(*M) = M_elements[1];
	int sum=0;
	for (int i=0; i<(Nd+2); i++){
		sum += N_elements[i];
	}
	(*N) = sum;
	double S_elements[Nd+2];
	S_elements[0] = 0.0;
	S_elements[1] = Sr;
	for (int i=2; i<(Nd+2); i++){
		S_elements[i] = Sd;
	}
	double L_elements[Nd+2];
	L_elements[0] = Lr;
	L_elements[1] = Lf;
	for (int i=2; i<(Nd+2); i++){
		L_elements[i] = Ld;
	}
	// Print antenna details
	printf(" n   N    M     L        S\n");
	printf("============================\n");
	for (int i=0; i<(Nd+2); i++){
		printf("%2d  %2d  %3d   %1.4f   %1.4f\n", i+1, N_elements[i], M_elements[i]+1, 
		L_elements[i], S_elements[i]);
	}
	printf("============================\n");
	printf("N = %d\n", (*N));
	
	double dL, L, S;
	segment* seg_ptr=malloc((*N)*sizeof(segment));
	
	int count=0;
	S=0.0;
	for (int i=0; i<(Nd+2); i++){
		L = L_elements[i];
		dL = L/(N_elements[i]+1.0);
		S += S_elements[i];
		for (int n=0; n<N_elements[i]; n++){
			seg_ptr[count].Index = i;
			// Points
			seg_ptr[count].rm.x = S;
			seg_ptr[count].rm.y = 0.0;
			seg_ptr[count].rm.z = (n+0)*dL-L/2.0;
			seg_ptr[count].rn.x = S;
			seg_ptr[count].rn.y = 0.0;
			seg_ptr[count].rn.z = (n+1)*dL-L/2.0;
			seg_ptr[count].rp.x = S;
			seg_ptr[count].rp.y = 0.0;
			seg_ptr[count].rp.z = (n+2)*dL-L/2.0;
			// Vectors
			seg_ptr[count].Lm.x = seg_ptr[count].rn.x-seg_ptr[count].rm.x;
			seg_ptr[count].Lm.y = seg_ptr[count].rn.y-seg_ptr[count].rm.y;
			seg_ptr[count].Lm.z = seg_ptr[count].rn.z-seg_ptr[count].rm.z;
			seg_ptr[count].Lm_hat = unitVector(seg_ptr[count].Lm);
			
			seg_ptr[count].Lp.x = seg_ptr[count].rp.x-seg_ptr[count].rn.x;
			seg_ptr[count].Lp.y = seg_ptr[count].rp.y-seg_ptr[count].rn.y;
			seg_ptr[count].Lp.z = seg_ptr[count].rp.z-seg_ptr[count].rn.z;
			seg_ptr[count].Lp_hat = unitVector(seg_ptr[count].Lp);
			count++;
		}
	}
	
	return seg_ptr;

}