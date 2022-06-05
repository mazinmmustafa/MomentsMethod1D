#include <stdlib.h>
// My header files
#include "Segment.h"

segment* VerticalDipole(int N, double L){
	segment* seg_ptr=malloc(N*sizeof(segment));
	double dL=L/(N+1.0);
	for (int n=0; n<N; n++){
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

segment* LShapeWire(int N1, int N2, double L1, double L2){
	segment* seg_ptr=malloc((N1+N2+1)*sizeof(segment));
	double dL1=L1/(N1+1.0);
	double dL2=L2/(N2+1.0);
	for (int n=0; n<N1; n++){
		// Points
		seg_ptr[n].rm.x = 0.0;
		seg_ptr[n].rm.y = (N1-n+1)*dL1;
		seg_ptr[n].rm.z = 0.0;
		seg_ptr[n].rn.x = 0.0;
		seg_ptr[n].rn.y = (N1-n+0)*dL1;
		seg_ptr[n].rn.z = 0.0;
		seg_ptr[n].rp.x = 0.0;
		seg_ptr[n].rp.y = (N1-n-1)*dL1;
		seg_ptr[n].rp.z = 0.0;
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
		// Corner segment
		// Points
		seg_ptr[N1].rm.x = 0.0;
		seg_ptr[N1].rm.y = +dL1;
		seg_ptr[N1].rm.z = 0.0;
		seg_ptr[N1].rn.x = 0.0;
		seg_ptr[N1].rn.y = 0.0;
		seg_ptr[N1].rn.z = 0.0;
		seg_ptr[N1].rp.x = 0.0;
		seg_ptr[N1].rp.y = 0.0;
		seg_ptr[N1].rp.z = +dL2;
		// Vectors
		seg_ptr[N1].Lm.x = seg_ptr[N1].rn.x-seg_ptr[N1].rm.x;
		seg_ptr[N1].Lm.y = seg_ptr[N1].rn.y-seg_ptr[N1].rm.y;
		seg_ptr[N1].Lm.z = seg_ptr[N1].rn.z-seg_ptr[N1].rm.z;
		seg_ptr[N1].Lm_hat = unitVector(seg_ptr[N1].Lm);
		
		seg_ptr[N1].Lp.x = seg_ptr[N1].rp.x-seg_ptr[N1].rn.x;
		seg_ptr[N1].Lp.y = seg_ptr[N1].rp.y-seg_ptr[N1].rn.y;
		seg_ptr[N1].Lp.z = seg_ptr[N1].rp.z-seg_ptr[N1].rn.z;
		seg_ptr[N1].Lp_hat = unitVector(seg_ptr[N1].Lp);
	for (int n=0; n<N2; n++){
		// Points
		seg_ptr[n+N1+1].rm.x = 0.0;
		seg_ptr[n+N1+1].rm.y = 0.0;
		seg_ptr[n+N1+1].rm.z = (n+0)*dL2;
		seg_ptr[n+N1+1].rn.x = 0.0;
		seg_ptr[n+N1+1].rn.y = 0.0;
		seg_ptr[n+N1+1].rn.z = (n+1)*dL2;
		seg_ptr[n+N1+11].rp.x = 0.0;
		seg_ptr[n+N1+1].rp.y = 0.0;
		seg_ptr[n+N1+1].rp.z = (n+2)*dL2;
		// Vectors
		seg_ptr[n+N1+1].Lm.x = seg_ptr[n+N1+1].rn.x-seg_ptr[n+N1+1].rm.x;
		seg_ptr[n+N1+1].Lm.y = seg_ptr[n+N1+1].rn.y-seg_ptr[n+N1+1].rm.y;
		seg_ptr[n+N1+1].Lm.z = seg_ptr[n+N1+1].rn.z-seg_ptr[n+N1+1].rm.z;
		seg_ptr[n+N1+1].Lm_hat = unitVector(seg_ptr[n+N1+1].Lm);
		
		seg_ptr[n+N1+1].Lp.x = seg_ptr[n+N1+1].rp.x-seg_ptr[n+N1+1].rn.x;
		seg_ptr[n+N1+1].Lp.y = seg_ptr[n+N1+1].rp.y-seg_ptr[n+N1+1].rn.y;
		seg_ptr[n+N1+1].Lp.z = seg_ptr[n+N1+1].rp.z-seg_ptr[n+N1+1].rn.z;
		seg_ptr[n+N1+1].Lp_hat = unitVector(seg_ptr[n+N1+1].Lp);
	}
	return seg_ptr;
}