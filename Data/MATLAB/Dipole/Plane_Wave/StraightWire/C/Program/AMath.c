#include <math.h>

double sinc(double x){
	if (x==0.0){
		return 1.0;
	}
	else{
		return sin(x)/x;
	}
}

double deg2rad(double theta){
	double pi=M_PI;
	return theta*pi/180.0;
};

double rad2deg(double theta){
	double pi=M_PI;
	return theta*180.0/pi;
};