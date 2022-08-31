#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

double Clebsch_(double j1,double m1,double j2,double m2,double j,double m);
double TLM_(double l,double m,double j1,double m1,double j2,double m2);
double unitTLM_(double l,double m,double j1,double m1,double j2,double m2);
int numberOfStates_(int spinCount, int *allSpinsTimesTwo);
void getIx_(double complex **operator, int spinIndex, int *allSpinsTimesTwo, int spinCount);
void getIy_(double complex **operator, int spinIndex, int *allSpinsTimesTwo, int spinCount);
void getIz_(double complex **operator, int spinIndex, int *allSpinsTimesTwo, int spinCount);
void getIp_(double complex **operator, int spinIndex, int *allSpinsTimesTwo, int spinCount);
void getIm_(double complex **operator, int spinIndex, int *allSpinsTimesTwo, int spinCount);
void getTlm_(double complex **operator, int spinIndex, int *allSpinsTimesTwo, int spinCount, int L, int M);
