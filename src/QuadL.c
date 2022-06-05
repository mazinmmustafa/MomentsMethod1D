#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <complex.h>
//
#include "../include/QuadL.h"
#include "../include/GaussLegendreRule.h"

int NQuad=0;
double *xL=NULL, *wL=NULL;

void setQuadL(int N){
    assert(N>0);
    assert(xL==NULL);
    assert(wL==NULL);
    NQuad = N;
    xL = malloc(N*sizeof(double));
    wL = malloc(N*sizeof(double));
    GaussLegendreRule(NQuad, xL, wL);
}

void unsetQuadL(){
    assert(xL!=NULL);
    assert(wL!=NULL);
    free(xL);
    free(wL);
    xL = NULL;
    wL = NULL;
	NQuad = 0;
}

complex double QuadL_1D(complex double func(double, void*),
        void *funcArgs, double a, double b){
    // \int_{a}^{b} func(x) dx
    assert(xL!=NULL);
    assert(wL!=NULL);
    assert(NQuad>0);
    //
	double hp=(b+a)/2.0;
	double hm=(b-a)/2.0;
	complex double sum=0.0;
	for (int i=0; i<NQuad; i++){
		sum+=hm*wL[i]*func(hm*xL[i]+hp, funcArgs);
	}
	return sum;
}

complex double QuadL_2D(complex double func(double, double, void*),
		void *funcArgs, double a2, double b2, double a1, double b1){
	// \int_{a2}^{b2} \int_{a1}^{b1} func(x, y) dx dy
	double hp1=(b1+a1)/2.0;
	double hm1=(b1-a1)/2.0;
	double hp2=(b2+a2)/2.0;
	double hm2=(b2-a2)/2.0;
	complex double sum=0.0;
	for (int i=0; i<NQuad; i++){
		for (int j=0; j<NQuad; j++){
			sum+=hm1*hm2*wL[i]*wL[j]*func(hm1*xL[i]+hp1, hm2*xL[j]+hp2, funcArgs);
		}
	}
	return sum;
}

complex double QuadL_Tri_2D(complex double func(double, double, void*), void *funcArgs){
	// \int_{0}^{1} \int_{0}^{1-\alpha} func(\alpha, \beta) d\beta d\alpha
	double zeta, eta;
	double alpha, beta, w;
	complex double sum=0.0;
	for (int i=0; i<NQuad; i++){
		zeta = xL[i];
		for (int j=0; j<NQuad; j++){
			eta = xL[j];
			alpha = (3.0+3.0*zeta-eta-zeta*eta)/8.0;
			beta = (3.0+3.0*eta-zeta-zeta*eta)/8.0;
			w = wL[i]*wL[j]*(2.0-zeta-eta)/16.0;
			sum+=w*func(alpha, beta, funcArgs);
		}
	}
	return sum;
}

complex double QuadL_Tri_4D(complex double func(double, double, double, double, void*),
        void *funcArgs){
	// \int_{0}^{1} \int_{0}^{1-\alpha2} \int_{0}^{1} \int_{0}^{1-\alpha1}
	// func(\alpha1, \beta1, \alpha2, \beta2) d\beta1 d\alpha1 d\beta2 d\alpha2
	double alpha1, beta1, zeta1, eta1, w1;
	double alpha2, beta2, zeta2, eta2, w2;
	complex double sum=0.0;
	for (int i=0; i<NQuad; i++){
		zeta2 = xL[i];
		for (int j=0; j<NQuad; j++){
			eta2 = xL[j];
			alpha2 = (3.0+3.0*zeta2-eta2-zeta2*eta2)/8.0;
			beta2 = (3.0+3.0*eta2-zeta2-zeta2*eta2)/8.0;
			w2 = wL[i]*wL[j]*(2.0-zeta2-eta2)/16.0;
			for (int ii=0; ii<NQuad; ii++){
				zeta1 = xL[ii];
				for (int jj=0; jj<NQuad; jj++){
					eta1 = xL[jj];
					alpha1 = (3.0+3.0*zeta1-eta1-zeta1*eta1)/8.0;
					beta1 = (3.0+3.0*eta1-zeta1-zeta1*eta1)/8.0;
					w1 = wL[ii]*wL[jj]*(2.0-zeta1-eta1)/16.0;
					sum+=w1*w2*func(alpha1, beta1, alpha2, beta2, funcArgs);
				}
			}
		}
	}
	return sum;
}
