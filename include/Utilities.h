#ifndef UTILITIES_H
#define UTILITIES_H

//
#include <stdio.h>
#include <time.h>

// Definitions
typedef struct Timer{
  time_t start, stop;
  double duration;
}Timer;
#define DefaultTimer {clock(), clock(), 0}

typedef struct LogFile{
  char *fileName;
  FILE *file;
}LogFile;
#define DefaultLogFile {NULL, NULL}

// Functions
void ticTimer(Timer *T);
void tocTimer(Timer *T);
void showTimer(Timer T);
void progressBar(int n, int N);
void openLogFile(LogFile *myLogFile, char *fileName);
void closeLogFile(LogFile *myLogFile);
void appendLogFile(LogFile myLogFile, char *newLine);
double deg2rad(double theta);
double rad2deg(double theta);
double roundn(double x, int n);

#endif