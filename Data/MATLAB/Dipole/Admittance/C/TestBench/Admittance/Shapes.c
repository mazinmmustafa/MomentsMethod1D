#include <stdlib.h>
// My header files
#include "Segment.h"

segment* VerticalDipoleDeltaGap(int* N, double L, int* M){
	if ((*N)%2==0){(*N)++;}
	(*M)=((*N)-1)/2;
	segment* seg_ptr=malloc((*N)*sizeof(segment));
	double dL=L/((*N)+1.0);
	for (int n=0; n<(*N); n++){
		// Points
		seg_ptr[n].rm.x = 0.0;
		seg_ptr[n].rm.y = 0.0;
		seg_ptr[n].rm.z = (n+0)*dL-L/2.0;
		seg_ptr[n].rn.x = 0.0;
		seg_ptr[n].rn.y = 0.0;
		seg_ptr[n].rn.z = (n+1)*dL-L/2.0;
		seg_ptr[n].rp.x = 0.0;
		seg_ptr[n].rp.y = 0.0;
		seg_ptr[n].rp.z = (n+2)*dL-L/2.0;
		// Vectors
		seg_ptr[n].Lm.x = seg_ptr[n].rn.x-seg_ptr[n].rm.x;
		seg_ptr[n].Lm.y = seg_ptr[n].rn.y-seg_ptr[n].rm.y;
		seg_ptr[n].Lm.z = seg_ptr[n].rn.z-seg_ptr[n].rm.z;
		seg_ptr[n].Lm_hat = unitVector(seg_ptr[n].Lm);
		
		seg_ptr[n].Lp.x = seg_ptr[n].rp.x-seg_ptr[n].rn.x;
		seg_ptr[n].Lp.y = seg_ptr[n].rp.y-seg_ptr[n].rn.y;
		seg_ptr[n].Lp.z = seg_ptr[n].rp.z-seg_ptr[n].rn.z;
		seg_ptr[n].Lp_hat = unitVector(seg_ptr[n].Lp);
	}
	return seg_ptr;
}