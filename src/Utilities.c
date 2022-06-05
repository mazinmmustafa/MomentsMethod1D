#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
//
#include "../include/Utilities.h"

void ticTimer(Timer *T){
  (*T).start = clock();
}

void tocTimer(Timer *T){
  (*T).stop = clock();
  (*T).duration = (double)((*T).stop-(*T).start)/CLOCKS_PER_SEC;
}

void showTimer(Timer T){
  double dT=T.duration/60.0; // Time in minutes
  if (dT<1.0){
    printf("Elapsed time is %0.2f seconds\n", dT*60.0);
  }else
  if (dT>=1.0&&dT<60.0){
    printf("Elapsed time is %0.2f minutes\n", dT);
  }else{
    printf("Elapsed time is %0.2f hours\n", dT/60.0);
  }
}

void progressBar(int n, int N){
  int NBar=30;
  n++;
  printf("\rProgess: [");
  for (int i=0; i<NBar*n/N; i++){
    printf("#");
  }
  for (int i=NBar*n/N; i<NBar; i++){
    printf("-");
  }
  printf("] %3d%%", 100*n/N);
  if (n==N){
    printf(" Done!\n");
  }
}

void openLogFile(LogFile *myLogFile, char *fileName){
  assert(fileName!=NULL);
  assert((*myLogFile).file==NULL);
  (*myLogFile).fileName = fileName;
  (*myLogFile).file = fopen((*myLogFile).fileName, "w");
  if ((*myLogFile).file==NULL){
    printf("Error: Unable to open '%s' file!", (*myLogFile).fileName);
    exit(1);
  }else{
    time_t timeStamp;
    time(&timeStamp);
    fprintf((*myLogFile).file, "%s\n", ctime(&timeStamp));
  }
}

void closeLogFile(LogFile *myLogFile){
  assert((*myLogFile).fileName!=NULL);
  assert((*myLogFile).file!=NULL);
  fclose((*myLogFile).file);
  (*myLogFile).fileName=NULL;
  (*myLogFile).file=NULL;
}

void appendLogFile(LogFile myLogFile, char *newLine){
  assert(myLogFile.fileName!=NULL);
  assert(myLogFile.file!=NULL);
  assert(newLine!=NULL);
  fprintf(myLogFile.file, "%s", newLine);
}

double deg2rad(double theta){
  return theta*(M_PI/180.0);
}

double rad2deg(double theta){
  return theta*(180.0/M_PI);
}

double roundn(double x, int n){
	double r=pow(10.0, n);
	return (double) round(x*r)/r;
}