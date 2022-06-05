#include <math.h>
#include <complex.h>
// My header files
#include "Type.h"
#include "Segment.h"

void FarField(int N, const segment* shape, double theta, double phi, 
           TYPE* In, double* E_theta_out, double* E_phi_out){
	TYPE j=I;
	double pi=M_PI;
	double k=2.0*pi;
	TYPE kn_p, kn_m, dot_p, dot_m;
	vector r_hat;
	r_hat.x = sin(theta)*cos(phi);
	r_hat.y = sin(theta)*sin(phi);
	r_hat.z = cos(theta);
	vector theta_hat, phi_hat;
	theta_hat.x =  cos(theta)*cos(phi);
	theta_hat.y =  cos(theta)*sin(phi);
	theta_hat.z = -sin(theta);
	phi_hat.x = -sin(phi);
	phi_hat.y =  cos(phi);
	phi_hat.z =  0.0;
	vector rn_m, rn_p;
	TYPE E_theta=(TYPE) 0.0, E_phi=(TYPE) 0.0;
	for (int n=0; n<N; n++){
		rn_m.x = shape[n].rm.x;
		rn_m.y = shape[n].rm.y;
		rn_m.z = shape[n].rm.z;
		rn_p.x = shape[n].rp.x;
		rn_p.y = shape[n].rp.y;
		rn_p.z = shape[n].rp.z;
		// kn+
		dot_p = dotVector(shape[n].Lp, r_hat);
		if (round(cabs(dot_p)*1E10)/1E10==0.0){
			kn_p = cexp(j*k*dotVector(rn_p, r_hat))/2.0;
		}
		else{
			kn_p = (cexp(-j*k*dot_p)*(1.0+j*k*dot_p)-1.0)/(k*k*dot_p*dot_p);
			kn_p *= cexp(j*k*dotVector(rn_p, r_hat));
		}
		// kn-
		dot_m = dotVector(shape[n].Lm, r_hat);
		if (round(cabs(dot_m)*1E10)/1E10==0.0){
			kn_m = cexp(j*k*dotVector(rn_m, r_hat))/2.0;
		}
		else{
			kn_m = (cexp(+j*k*dot_m)*(1.0-j*k*dot_m)-1.0)/(k*k*dot_m*dot_m);
			kn_m *= cexp(j*k*dotVector(rn_m, r_hat));
		}
		E_theta += In[n]*(kn_p*dotVector(theta_hat, shape[n].Lp)+
		                      kn_m*dotVector(theta_hat, shape[n].Lm));
		E_phi += In[n]*(kn_p*dotVector(phi_hat, shape[n].Lp)+
		                    kn_m*dotVector(phi_hat, shape[n].Lm));
	}
	E_theta = cabs(E_theta);
	E_phi = cabs(E_phi);
	*E_theta_out = E_theta;
	*E_phi_out = E_phi;
}