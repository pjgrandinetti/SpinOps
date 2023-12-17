#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

typedef struct EULER {
	double alpha;				/* Euler angle alpha */
	double beta;				/* Euler angle beta */
	double gamma;				/* Euler angle gamma */
	} euler;

void getrho2_pas_(double complex *tensor, double zeta, double eta);
void getrho1_pas_(double complex *tensor, double zeta);
double wigner_d_(double l,double m1,double m2,double beta);
double complex DLM_(double l,double  m1,double m2, double alpha, double beta, double gamma);
void Rot_(double j, double complex *initial, double alpha, double beta, double gamma, double complex *final);
