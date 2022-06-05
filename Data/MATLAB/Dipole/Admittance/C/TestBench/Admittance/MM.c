#include <math.h>
// My header files
#include "Type.h"
#include "Segment.h"
#include "Quad.h"
#include "Integrands.h"
#include "Matrix.h"
#include "Singularity.h"

void MM_DeltaGap(int N, int M, double a, TYPE** Z, TYPE* V, const segment* shape){
	TYPE j=I;
	double pi=M_PI;
	double eta=120*pi;
	double k=2.0*pi;
	TYPE term1, term2;
	for (int m=0; m<N; m++){
		for (int n=0; n<N; n++){
			if (m==n){	
				term1 =  I1(magVector(shape[m].Lp), k, a) 
				        +psi_pm(shape[m], shape[n], k, a)
						+psi_mp(shape[m], shape[n], k, a)
						+I1(magVector(shape[m].Lm), k, a);
				term2 =  I2(magVector(shape[m].Lp), k, a) 
				        -phi_pm(shape[m], shape[n], k, a)
						-phi_mp(shape[m], shape[n], k, a)
						+I2(magVector(shape[m].Lm), k, a);
			}
			else if (m==n+1){
				term1 =  psi_pp(shape[m], shape[n], k, a)
				        +psi_pm(shape[m], shape[n], k, a)
						+I3(magVector(shape[m].Lm), k, a)
						+psi_mm(shape[m], shape[n], k, a);
				term2 =  phi_pp(shape[m], shape[n], k, a)
				        -phi_pm(shape[m], shape[n], k, a)
						-I2(magVector(shape[m].Lm), k, a)
						+phi_mm(shape[m], shape[n], k, a);
			}
			else if (n==m+1){
				term1 =  psi_pp(shape[m], shape[n], k, a)
				        +I3(magVector(shape[m].Lp), k, a)
				        +psi_mp(shape[m], shape[n], k, a)
						+psi_mm(shape[m], shape[n], k, a);
				term2 =  phi_pp(shape[m], shape[n], k, a)
				        -I2(magVector(shape[m].Lp), k, a)
				        -phi_mp(shape[m], shape[n], k, a)
						+phi_mm(shape[m], shape[n], k, a);
			}
			else{
				term1 =  psi_pp(shape[m], shape[n], k, a)
				        +psi_pm(shape[m], shape[n], k, a)
				        +psi_mp(shape[m], shape[n], k, a)
						+psi_mm(shape[m], shape[n], k, a);
				term2 =  phi_pp(shape[m], shape[n], k, a)
				        -phi_pm(shape[m], shape[n], k, a)
				        -phi_mp(shape[m], shape[n], k, a)
						+phi_mm(shape[m], shape[n], k, a);
			}
			Z[m][n] = j*k*eta*term1-(j*eta/k)*term2;
		}
		if (m==M){
			V[m] = 1.0;
		}
	}
}