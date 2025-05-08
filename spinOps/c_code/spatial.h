#ifndef SPATIAL_H
#define SPATIAL_H

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h> 
#include <float.h>
#include <stdbool.h>

typedef struct EULER {
	double alpha;				/* Euler angle alpha */
	double beta;				/* Euler angle beta */
	double gamma;				/* Euler angle gamma */
	} euler;

/*!
 @function wigner_d_
 @abstract Computes the Wigner small-d function d(l, m1, m2, beta).
 @discussion The Wigner small-d function is a component of the Wigner rotation matrix element 
              and depends on the rank `l`, the orders 
              `m1` and `m2`, and the Euler angle `beta`. For `l = 2`, the function uses explicit 
              formulas for efficiency. For general `l`, it computes the value using a summation 
              formula involving factorials and powers of trigonometric functions.
 @param two_l 2 times the rank (non-negative integer or half-integer).
 @param two_m1 2 times the order in the initial frame (-l <= m1 <= l).
 @param two_m2 2 times the order in the final frame (-l <= m2 <= l).
 @param beta The second Euler angle (rotation about the y-axis).
 @return The value of the Wigner small-d function d(l, m1, m2, beta).
 */
double wigner_d_(const int two_J, const int two_m1, const int two_m2, const double theta);

/*!
 @function DLM_
 @abstract Computes the Wigner rotation matrix element D(l, m1, m2).
 @discussion The Wigner rotation matrix element is defined as:
              D(l, m1, m2) = e^(-i * m1 * alpha) * d(l, m1, m2, beta) * e^(-i * m2 * gamma),
              where d(l, m1, m2, beta) is the Wigner small-d function, and
              alpha, beta, gamma are the Euler angles.
 @param two_l Two times the rank (non-negative integer or half-integer).
 @param two_m1 Two times the order in the initial frame (-l <= m1 <= l).
 @param two_m2 Two times the order in the final frame (-l <= m2 <= l).
 @param alpha The first Euler angle (rotation about the z-axis).
 @param beta The second Euler angle (rotation about the y-axis).
 @param gamma The third Euler angle (rotation about the z-axis).
 @return The complex value of the Wigner rotation matrix element D(l, m1, m2).
 */
double complex DLM_(const int two_l, const int two_m1, const int two_m2, const double alpha, const double beta, const double gamma);
/*!
 @function Rot_
 @abstract Performs a rotational transformation of a spherical tensor from one frame to another.
 @discussion This function applies a rotational transformation to a tensor represented 
              in the initial frame, using the Euler angles `alpha`, `beta`, and `gamma`. 
              The transformation is performed using Wigner rotation matrix elements. 
              For `j = 2`, the function uses the `wigner_d_` function for efficiency. 
              For general `j`, it uses the `DLM_` function to compute the rotation matrix.
 @param two_j The spherical tensor of the tensor (non-negative integer or half-integer) times 2.
 @param initial A pointer to the array representing the tensor components in the initial frame.
 @param alpha The first Euler angle (rotation about the z-axis).
 @param beta The second Euler angle (rotation about the y-axis).
 @param gamma The third Euler angle (rotation about the z-axis).
 @param final A pointer to the array where the transformed spherical tensor components will be stored.
 */
void Rot_(const int two_j, const double complex *initial,const double alpha, const double beta, const double gamma,double complex *final);


/*!
 @function get_rho2_pas_
 @abstract Creates the complex vector for the spherical tensor of rank 2.
 @discussion The spherical tensor of rank 2 has five components corresponding to 
              the orders m = -2, -1, 0, +1, +2. These components 
              are mapped to the indices 0, 1, 2, 3, and 4 in the array. The function 
              initializes the tensor based on the traceless 2nd-rank symmetric tensor 
              anisotropy `zeta` and the asymmetry parameter `eta`.
 @param tensor A pointer to the array representing the spherical tensor components.
 @param zeta The traceless 2nd-rank symmetric tensor anisotropy.
 @param eta The traceless 2nd-rank symmetric tensor asymmetry parameter.
 */
void get_rho2_pas_(double complex *tensor, const double zeta, const double eta);
/*!
 @function get_rho1_pas_
 @abstract Creates the complex vector for the spherical tensor of rank 1.
 @discussion The spherical tensor of rank 1 has three components corresponding to 
              the orders m = -1, 0, +1. These components are mapped 
              to the indices 0, 1, and 2 in the array. The function initializes the 
              tensor based on the traceless 1st-rank symmetric tensor anisotropy `zeta`.
 @param tensor A pointer to the array representing the spherical tensor components.
 @param zeta The traceless 1st-rank symmetric tensor anisotropy.
 */
void get_rho1_pas_(double complex *tensor, const double zeta);



#endif // SPATIAL_H
