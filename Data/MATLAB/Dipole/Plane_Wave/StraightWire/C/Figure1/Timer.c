#include <stdio.h>
#include <time.h>

clock_t Start, Stop;

void tic(){
	Start = clock();
}

void toc(){
	Stop = clock();
	double time=(double) (Stop-Start)/CLOCKS_PER_SEC;
	printf("\n========================================\n");
	if (time>=60&&time<3600){
		printf("Elapsed time %1.2f minutes\n", time/60);
	}
	else if (time>=3600){
		printf("Elapsed time %1.2f hours\n", time/3600);
	}
	else{
		printf("Elapsed time %1.2f seconds\n", time);
	}
	printf("========================================\n");
}