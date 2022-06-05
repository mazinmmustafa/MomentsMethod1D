#include <stdio.h> 
#include <math.h>
#include <complex.h>
// My header files
#include "Type.h"
#include "Singularity.h"

int main(){
	
	double pi=M_PI;
	double a, L, k;
	
	k = 2.0*pi;
	L = 0.05;
	a = 1.0E-4;
	
	TYPE ans;
	ans = I1(L, k, a);
	printf("I1 = (%12.4E, %12.4E)\n", creal(ans), cimag(ans));
	ans = I2(L, k, a);
	printf("I2 = (%12.4E, %12.4E)\n", creal(ans), cimag(ans));
	ans = I3(L, k, a);
	printf("I3 = (%12.4E, %12.4E)\n", creal(ans), cimag(ans));
	
}

