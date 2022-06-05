#include <math.h>
#include <complex.h>
// My header files
#include "Type.h"
#include "Segment.h"

void sigma(int N, const segment* shape, double theta, double phi, 
           TYPE* In, double* sigma_theta_out, double* sigma_phi_out){
	TYPE j=I;
	double pi=M_PI;
	double k=2.0*pi;
	double eta=120.0*pi;
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
	TYPE sigma_theta=(TYPE) 0.0, sigma_phi=(TYPE) 0.0;
	for (int n=0; n<N; n++){
		rn_m.x = shape[n].rm.x;
		rn_m.y = shape[n].rm.y;
		rn_m.z = shape[n].rm.z;
		rn_p.x = shape[n].rp.x;
		rn_p.y = shape[n].rp.y;
		rn_p.z = shape[n].rp.z;
		// kn+
		dot_p = dotVector(shape[n].Lp, r_hat);
		if (cabs(dot_p)==0.0){
			kn_p = cexp(j*k*dotVector(rn_p, r_hat))/2.0;
		}
		else{
			kn_p = (cexp(-j*k*dot_p)*(1.0+j*k*dot_p)-1.0)/(k*k*dot_p*dot_p);
			kn_p *= cexp(j*k*dotVector(rn_p, r_hat));
		}
		// kn-
		dot_m = dotVector(shape[n].Lm, r_hat);
		if (cabs(dot_m)==0.0){
			kn_m = cexp(j*k*dotVector(rn_m, r_hat))/2.0;
		}
		else{
			kn_m = (cexp(+j*k*dot_m)*(1.0-j*k*dot_m)-1.0)/(k*k*dot_m*dot_m);
			kn_m *= cexp(j*k*dotVector(rn_m, r_hat));
		}
		sigma_theta += In[n]*(kn_p*dotVector(theta_hat, shape[n].Lp)+
		                      kn_m*dotVector(theta_hat, shape[n].Lm));
		sigma_phi += In[n]*(kn_p*dotVector(phi_hat, shape[n].Lp)+
		                    kn_m*dotVector(phi_hat, shape[n].Lm));
	}
	sigma_theta = pi*eta*eta*cabs(sigma_theta)*cabs(sigma_theta);
	sigma_phi = pi*eta*eta*cabs(sigma_phi)*cabs(sigma_phi);
	*sigma_theta_out = sigma_theta;
	*sigma_phi_out = sigma_phi;
}