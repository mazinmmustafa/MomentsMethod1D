#include <complex.h>
#include <math.h>
// My header files
#include "Type.h"
#include "Segment.h"
#include "Quad.h"

TYPE psi_pp(segment seg_m, segment seg_n, double k, double a){
	TYPE j=I;
	double pi=M_PI;
	double R(TYPE alpha, TYPE alpha_){
		double x, y, z;
		x = cabs(seg_m.rp.x-seg_n.rp.x-alpha*seg_m.Lp.x+alpha_*seg_n.Lp.x);
		y = cabs(seg_m.rp.y-seg_n.rp.y-alpha*seg_m.Lp.y+alpha_*seg_n.Lp.y);
		z = cabs(seg_m.rp.z-seg_n.rp.z-alpha*seg_m.Lp.z+alpha_*seg_n.Lp.z);
		return sqrt(x*x+y*y+z*z+a*a);
	};
	TYPE g(TYPE alpha, TYPE alpha_){
		return cexp(-j*k*R(alpha, alpha_))/(4.0*pi*R(alpha, alpha_));
	};
	TYPE integrand(TYPE alpha, TYPE alpha_){
		return alpha*alpha_*g(alpha, alpha_);
	};
	return dotVector(seg_m.Lp, seg_n.Lp)*Quad32_2D(integrand, 0.0, 1.0, 0.0, 1.0);
}

TYPE psi_pm(segment seg_m, segment seg_n, double k, double a){
	TYPE j=I;
	double pi=M_PI;
	double R(TYPE alpha, TYPE alpha_){
		double x, y, z;
		x = cabs(seg_m.rp.x-seg_n.rm.x-alpha*seg_m.Lp.x-alpha_*seg_n.Lm.x);
		y = cabs(seg_m.rp.y-seg_n.rm.y-alpha*seg_m.Lp.y-alpha_*seg_n.Lm.y);
		z = cabs(seg_m.rp.z-seg_n.rm.z-alpha*seg_m.Lp.z-alpha_*seg_n.Lm.z);
		return sqrt(x*x+y*y+z*z+a*a);
	};
	TYPE g(TYPE alpha, TYPE alpha_){
		return cexp(-j*k*R(alpha, alpha_))/(4.0*pi*R(alpha, alpha_));
	};
	TYPE integrand(TYPE alpha, TYPE alpha_){
		return alpha*alpha_*g(alpha, alpha_);
	};
	return dotVector(seg_m.Lp, seg_n.Lm)*Quad32_2D(integrand, 0.0, 1.0, 0.0, 1.0);
}

TYPE psi_mp(segment seg_m, segment seg_n, double k, double a){
	TYPE j=I;
	double pi=M_PI;
	double R(TYPE alpha, TYPE alpha_){
		double x, y, z;
		x = cabs(seg_m.rm.x-seg_n.rp.x+alpha*seg_m.Lm.x+alpha_*seg_n.Lp.x);
		y = cabs(seg_m.rm.y-seg_n.rp.y+alpha*seg_m.Lm.y+alpha_*seg_n.Lp.y);
		z = cabs(seg_m.rm.z-seg_n.rp.z+alpha*seg_m.Lm.z+alpha_*seg_n.Lp.z);
		return sqrt(x*x+y*y+z*z+a*a);
	};
	TYPE g(TYPE alpha, TYPE alpha_){
		return cexp(-j*k*R(alpha, alpha_))/(4.0*pi*R(alpha, alpha_));
	};
	TYPE integrand(TYPE alpha, TYPE alpha_){
		return alpha*alpha_*g(alpha, alpha_);
	};
	return dotVector(seg_m.Lm, seg_n.Lp)*Quad32_2D(integrand, 0.0, 1.0, 0.0, 1.0);
}

TYPE psi_mm(segment seg_m, segment seg_n, double k, double a){
	TYPE j=I;
	double pi=M_PI;
	double R(TYPE alpha, TYPE alpha_){
		double x, y, z;
		x = cabs(seg_m.rm.x-seg_n.rm.x+alpha*seg_m.Lm.x-alpha_*seg_n.Lm.x);
		y = cabs(seg_m.rm.y-seg_n.rm.y+alpha*seg_m.Lm.y-alpha_*seg_n.Lm.y);
		z = cabs(seg_m.rm.z-seg_n.rm.z+alpha*seg_m.Lm.z-alpha_*seg_n.Lm.z);
		return sqrt(x*x+y*y+z*z+a*a);
	};
	TYPE g(TYPE alpha, TYPE alpha_){
		return cexp(-j*k*R(alpha, alpha_))/(4.0*pi*R(alpha, alpha_));
	};
	TYPE integrand(TYPE alpha, TYPE alpha_){
		return alpha*alpha_*g(alpha, alpha_);
	};
	return dotVector(seg_m.Lm, seg_n.Lm)*Quad32_2D(integrand, 0.0, 1.0, 0.0, 1.0);
}

TYPE phi_pp(segment seg_m, segment seg_n, double k, double a){
	TYPE j=I;
	double pi=M_PI;
	double R(TYPE alpha, TYPE alpha_){
		double x, y, z;
		x = cabs(seg_m.rp.x-seg_n.rp.x-alpha*seg_m.Lp.x+alpha_*seg_n.Lp.x);
		y = cabs(seg_m.rp.y-seg_n.rp.y-alpha*seg_m.Lp.y+alpha_*seg_n.Lp.y);
		z = cabs(seg_m.rp.z-seg_n.rp.z-alpha*seg_m.Lp.z+alpha_*seg_n.Lp.z);
		return sqrt(x*x+y*y+z*z+a*a);
	};
	TYPE g(TYPE alpha, TYPE alpha_){
		return cexp(-j*k*R(alpha, alpha_))/(4.0*pi*R(alpha, alpha_));
	};
	TYPE integrand(TYPE alpha, TYPE alpha_){
		return g(alpha, alpha_);
	};
	return Quad32_2D(integrand, 0.0, 1.0, 0.0, 1.0);
}

TYPE phi_pm(segment seg_m, segment seg_n, double k, double a){
	TYPE j=I;
	double pi=M_PI;
	double R(TYPE alpha, TYPE alpha_){
		double x, y, z;
		x = cabs(seg_m.rp.x-seg_n.rm.x-alpha*seg_m.Lp.x-alpha_*seg_n.Lm.x);
		y = cabs(seg_m.rp.y-seg_n.rm.y-alpha*seg_m.Lp.y-alpha_*seg_n.Lm.y);
		z = cabs(seg_m.rp.z-seg_n.rm.z-alpha*seg_m.Lp.z-alpha_*seg_n.Lm.z);
		return sqrt(x*x+y*y+z*z+a*a);
	};
	TYPE g(TYPE alpha, TYPE alpha_){
		return cexp(-j*k*R(alpha, alpha_))/(4.0*pi*R(alpha, alpha_));
	};
	TYPE integrand(TYPE alpha, TYPE alpha_){
		return g(alpha, alpha_);
	};
	return Quad32_2D(integrand, 0.0, 1.0, 0.0, 1.0);
}

TYPE phi_mp(segment seg_m, segment seg_n, double k, double a){
	TYPE j=I;
	double pi=M_PI;
	double R(TYPE alpha, TYPE alpha_){
		double x, y, z;
		x = cabs(seg_m.rm.x-seg_n.rp.x+alpha*seg_m.Lm.x+alpha_*seg_n.Lp.x);
		y = cabs(seg_m.rm.y-seg_n.rp.y+alpha*seg_m.Lm.y+alpha_*seg_n.Lp.y);
		z = cabs(seg_m.rm.z-seg_n.rp.z+alpha*seg_m.Lm.z+alpha_*seg_n.Lp.z);
		return sqrt(x*x+y*y+z*z+a*a);
	};
	TYPE g(TYPE alpha, TYPE alpha_){
		return cexp(-j*k*R(alpha, alpha_))/(4.0*pi*R(alpha, alpha_));
	};
	TYPE integrand(TYPE alpha, TYPE alpha_){
		return g(alpha, alpha_);
	};
	return Quad32_2D(integrand, 0.0, 1.0, 0.0, 1.0);
}

TYPE phi_mm(segment seg_m, segment seg_n, double k, double a){
	TYPE j=I;
	double pi=M_PI;
	double R(TYPE alpha, TYPE alpha_){
		double x, y, z;
		x = cabs(seg_m.rm.x-seg_n.rm.x+alpha*seg_m.Lm.x-alpha_*seg_n.Lm.x);
		y = cabs(seg_m.rm.y-seg_n.rm.y+alpha*seg_m.Lm.y-alpha_*seg_n.Lm.y);
		z = cabs(seg_m.rm.z-seg_n.rm.z+alpha*seg_m.Lm.z-alpha_*seg_n.Lm.z);
		return sqrt(x*x+y*y+z*z+a*a);
	};
	TYPE g(TYPE alpha, TYPE alpha_){
		return cexp(-j*k*R(alpha, alpha_))/(4.0*pi*R(alpha, alpha_));
	};
	TYPE integrand(TYPE alpha, TYPE alpha_){
		return g(alpha, alpha_);
	};
	return Quad32_2D(integrand, 0.0, 1.0, 0.0, 1.0);
}
