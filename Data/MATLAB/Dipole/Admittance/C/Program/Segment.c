#include <stdio.h>
#include <math.h>
// My header files
#include "Segment.h"

double dotVector(vector v1, vector v2){
	return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

double magVector(vector v){
	return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}

vector unitVector(vector v){
	vector v_hat;
	double v_mag=magVector(v);
	v_hat.x = v.x/v_mag;
	v_hat.y = v.y/v_mag;
	v_hat.z = v.z/v_mag;
	return v_hat;
}

void printShape(int N, segment* seg_ptr){
	printf("Shape data:\n");
	printf("  n   |");
	printf("         rm          |");
	printf("         rn          |");
	printf("         rp           ");
	printf("\n========================================================================\n");
	for (int n=0; n<N; n++){
		printf("%3d   |", n+1);
		printf("(%5.2f, %5.2f, %5.2f) ", seg_ptr[n].rm.x, seg_ptr[n].rm.y, seg_ptr[n].rm.z);
		printf("(%5.2f, %5.2f, %5.2f) ", seg_ptr[n].rn.x, seg_ptr[n].rn.y, seg_ptr[n].rn.z);
		printf("(%5.2f, %5.2f, %5.2f) ", seg_ptr[n].rp.x, seg_ptr[n].rp.y, seg_ptr[n].rp.z);
		printf("\n");
	}
}