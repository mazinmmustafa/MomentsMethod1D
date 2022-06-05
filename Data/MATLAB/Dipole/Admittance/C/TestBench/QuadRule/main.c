#include <stdio.h> 
#include <math.h>
#include <complex.h>
// My header files
#include "Type.h"
#include "Quad.h"

int main(){
	
	TYPE func_1(TYPE x){
		return ccos(10*x+I);
	}
	
	TYPE ans1=Quad32_1D(func_1, 0.0, 2.0);
	printf("I = %1.14E +j %1.14E\n", creal(ans1), cimag(ans1));
	
	TYPE func_2(TYPE x, TYPE y){
		return ccos(10*x+I*y);
	}
	
	TYPE ans2=Quad32_2D(func_2, 0.0, 2.0, 0.0, 3.0);
	printf("I = %1.14E +j %1.14E\n", creal(ans2), cimag(ans2));
	
}

