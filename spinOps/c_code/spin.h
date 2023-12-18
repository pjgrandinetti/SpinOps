#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

double clebsch_(double j1,double m1,double j2,double m2,double j,double m);
double tlm_(double l,double m,double j1,double m1,double j2,double m2);
double unit_tlm_(double l,double m,double j1,double m1,double j2,double m2);
int numberOfStates_(int spinCount, int *spinsTimesTwo);
void getIx_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount);
void getIy_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount);
void getIz_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount);
void getIp_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount);
void getIm_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount);
void getTlm_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount, int L, int M);
void getTlm_unit_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount, int L, int M);

void getEf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);
void getIxf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);
void getIyf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);
void getIzf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);
void getIpf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);
void getImf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);


double mypow(double x, int n);
double fac(double x);
