#include <math.h>

double sinc(double x){
	if (x==0.0){
		return 1.0;
	}
	else{
		return sin(x)/x;
	}
}