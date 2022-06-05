#include <stdio.h>
#include <math.h>
// My header files
#include "Type.h"
#include "Segment.h"
#include "Quad.h"
#include "Integrands.h"
#include "Matrix.h"
#include "Singularity.h"

void MM_PlaneWave(int N, double a, TYPE E_TM, TYPE E_TE, 
                  double theta_i, double phi_i,  
                   TYPE** Z, TYPE* V, const segment* shape){
	TYPE j=I;
	double pi=M_PI;
	double eta=120*pi;
	double k=2.0*pi;
	vector theta_i_hat, phi_i_hat, k_i;
	theta_i_hat.x =  cos(theta_i)*cos(phi_i);
	theta_i_hat.y =  cos(theta_i)*sin(phi_i);
	theta_i_hat.z = -sin(theta_i);
	phi_i_hat.x = -sin(phi_i);
	phi_i_hat.y =  cos(phi_i);
	phi_i_hat.z =  0.0;
	k_i.x = k*sin(theta_i)*cos(phi_i);
	k_i.y = k*sin(theta_i)*sin(phi_i);
	k_i.z = k*cos(theta_i);
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
		TYPE integrand_m(TYPE alpha){
			return alpha*cexp(+j*alpha*dotVector(k_i, shape[m].Lm));
		};
		TYPE integrand_p(TYPE alpha){
			return alpha*cexp(-j*alpha*dotVector(k_i, shape[m].Lp));
		};
		TYPE chi_m, chi_p;
		vector rm_m, rm_p; 
		rm_m.x = shape[m].rm.x;
		rm_m.y = shape[m].rm.y;
		rm_m.z = shape[m].rm.z;
		rm_p.x = shape[m].rp.x;
		rm_p.y = shape[m].rp.y;
		rm_p.z = shape[m].rp.z;
		chi_m = cexp(j*dotVector(k_i, rm_m))*Quad32_1D(integrand_m, 0.0, 1.0);
		chi_p = cexp(j*dotVector(k_i, rm_p))*Quad32_1D(integrand_p, 0.0, 1.0);
		V[m] = E_TM*(dotVector(shape[m].Lp, theta_i_hat)*chi_p+  
		             dotVector(shape[m].Lm, theta_i_hat)*chi_m)
			  +E_TE*(dotVector(shape[m].Lp, phi_i_hat)*chi_p+  
		             dotVector(shape[m].Lm, phi_i_hat)*chi_m);
	}
}