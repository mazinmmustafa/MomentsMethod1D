#ifndef SEGMENT_H
#define SEGMENT_H

typedef struct{
	double x, y, z;
}vector;

typedef struct{
	double x, y, z;
}point;

typedef struct{
	int Index;
	point rm, rn, rp;
	vector Lm, Lp, Lm_hat, Lp_hat;
}segment;

// Functions
double dotVector(vector, vector);
double magVector(vector);
vector unitVector(vector);
void printShape(int, segment*);

#endif