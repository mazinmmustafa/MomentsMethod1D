#include <math.h>
#include <complex.h>
// My header files
#include "Type.h"
#include "Quad.h"
#include "AMath.h"

TYPE I1(double L, double k, double a){
	TYPE j=I;
	double pi=M_PI;
	TYPE term1, term2, term3;
	double R(TYPE alpha, TYPE alpha_){
		return sqrt(L*L*(alpha-alpha_)*(alpha-alpha_)+a*a);
	};
	TYPE func1(TYPE alpha, TYPE alpha_){
		return alpha*alpha_*cexp(-j*k*R(alpha, alpha_)/2.0)*sinc(k*R(alpha, alpha_)/2.0);
	};
	TYPE func2(TYPE alpha){
		double A=sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a)+(1.0-alpha)*L; 
		double B=sqrt(alpha*alpha*L*L+a*a)-alpha*L;
		return alpha*alpha*log(A/B);
	};
	TYPE func3(TYPE alpha){
		double A=sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a); 
		double B=sqrt(alpha*alpha*L*L+a*a);
		return alpha*(A-B);
	};
	term1 = (-j*k*L*L/(4.0*pi))*Quad32_2D(func1, 0.0, 1.0, 0.0, 1.0);
	term2 = (L/(4.0*pi))*Quad32_1D(func2, 0.0, 1.0);
	term3 = (1.0/(4.0*pi))*Quad32_1D(func3, 0.0, 1.0);
	return term1+term2+term3;
}

TYPE I2(double L, double k, double a){
	TYPE j=I;
	double pi=M_PI;
	TYPE term1, term2;
	double R(TYPE alpha, TYPE alpha_){
		return sqrt(L*L*(alpha-alpha_)*(alpha-alpha_)+a*a);
	};
	TYPE func1(TYPE alpha, TYPE alpha_){
		return cexp(-j*k*R(alpha, alpha_)/2.0)*sinc(k*R(alpha, alpha_)/2.0);
	};
	TYPE func2(TYPE alpha){
		double A=sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a)+(1.0-alpha)*L; 
		double B=sqrt(alpha*alpha*L*L+a*a)-alpha*L;
		return log(A/B);
	};
	term1 = (-j*k/(4.0*pi))*Quad32_2D(func1, 0.0, 1.0, 0.0, 1.0);
	term2 = (1.0/(4.0*pi*L))*Quad32_1D(func2, 0.0, 1.0);
	return term1+term2;
}

TYPE I3(double L, double k, double a){
	TYPE j=I;
	double pi=M_PI;
	TYPE term1, term2, term3;
	double R(TYPE alpha, TYPE alpha_){
		return sqrt(L*L*(1.0-alpha-alpha_)*(1.0-alpha-alpha_)+a*a);
	};
	TYPE func1(TYPE alpha, TYPE alpha_){
		return alpha*alpha_*cexp(-j*k*R(alpha, alpha_)/2.0)*sinc(k*R(alpha, alpha_)/2.0);
	};
	TYPE func2(TYPE alpha){
		double A=sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a)+(1.0-alpha)*L; 
		double B=sqrt(alpha*alpha*L*L+a*a)-alpha*L;
		return alpha*(1.0-alpha)*log(A/B);
	};
	TYPE func3(TYPE alpha){
		double A=sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a); 
		double B=sqrt(alpha*alpha*L*L+a*a);
		return alpha*(A-B);
	};
	term1 = (-j*k*L*L/(4.0*pi))*Quad32_2D(func1, 0.0, 1.0, 0.0, 1.0);
	term2 = (L/(4.0*pi))*Quad32_1D(func2, 0.0, 1.0);
	term3 = (1.0/(4.0*pi))*Quad32_1D(func3, 0.0, 1.0);
	return term1+term2-term3;
}