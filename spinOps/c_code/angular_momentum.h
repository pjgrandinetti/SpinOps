#ifndef ANGULAR_MOMENTUM_H
#define ANGULAR_MOMENTUM_H

/* angular_momentum.h */
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h> 
#include <float.h>
#include <stdbool.h>

#define MAX_SMALL_FAC 32

// Maximum supported spin twice-I and tensor rank
#define MAX_TWO_I  11    // supports 2I = 1,2,…,11  i.e. I = ½,1,…,11/2
#define MAX_L       8    // supports l = 0,1,…,8

// Utility macros
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

extern const double small_fac[MAX_SMALL_FAC + 1];

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


/*!
 @function clebsch_
 @abstract Calculates the Clebsch-Gordon coefficients.
 @discussion This function computes the Clebsch-Gordon coefficients 
              `< j, m | j1, j2, m1, m2 >`. The Clebsch-Gordon coefficients 
              are used in quantum mechanics to describe the coupling of angular 
              momenta. The function ensures that the input values satisfy the 
              necessary conditions for valid coefficients.
 @param two_J1 The first angular momentum quantum number, j1 times 2.
 @param two_M1 The projection quantum number associated with `j1` times 2.
 @param two_J2 The second angular momentum quantum number times 2.
 @param two_M2 The projection quantum number associated with `j2` times 2.
 @param two_J The total angular momentum quantum number, `j` times 2.
 @param two_M The total projection quantum number times 2.
 @return The Clebsch-Gordon coefficient `< j, m | j1, j2, m1, m2 >` as a double. 
         Returns 0 if the input values do not satisfy the necessary conditions.
 */
double clebsch_(const int two_J1, const int two_M1,const int two_J2, const int two_M2,const int two_J,  const int two_M);

/*!
 @function tlm_
 @abstract Evaluates the matrix element `⟨I,m1| T_{l,m} |I,m2⟩`.
 @discussion This function calculates the matrix element `⟨I,m1| T_{l,m} |I,m2⟩` 
              using the definition from Bowden and Hutchinson, J. Magn. Reson. 67, 403, 1986. 
              The calculation involves Clebsch-Gordon coefficients and reduced matrix elements. 
 @param l       The rank of the tensor operator.
 @param m       The order of the tensor operator.
 @param two_I   2×I
 @param two_m1  The projection quantum number associated with I.
 @param two_m2  The projection quantum number associated with I.
 @return        The matrix element `⟨I,m1| T_{l,m} |I,m2⟩` as a double. 
 */
double tlm_(const int l,
            const int m,
            const int two_I,
            const int two_m1,
            const int two_m2);

/*!
 @function unit_tlm_
 @abstract Evaluates the matrix element `<j1 m1|T_hat(lm)|j2 m2>` for unit tensors.
 @discussion This function calculates the matrix element `<j1 m1|T_hat(lm)|j2 m2>` 
              using the definition of unit tensors from Bowden and Hutchinson, 
              J. Magn. Reson. 67, 403, 1986. The calculation involves Clebsch-Gordon 
              coefficients and normalization factors. The function assumes that 
              `j1` equals `j2` for the calculation.
 @param l The rank of the tensor operator.
 @param m The order of the tensor operator.
 @param two_j1 The first angular momentum quantum number times 2.
 @param two_m1 The projection quantum number associated with `j1` times 2.
 @param two_j2 The second angular momentum quantum number times 2.
 @param two_m2 The projection quantum number associated with `j2` times 2.
 @return The matrix element `<j1 m1|T_hat(lm)|j2 m2>` as a double. Returns 0 if `j1` is not equal to `j2`.
 */
double unit_tlm_(const int l,
                 const int m,
                 const int two_I,
                 const int two_m1,
                 const int two_m2);

/*!
 @function number_of_states_
 @abstract Calculates the size of the state space for a spin system.
 @discussion This function computes the total number of quantum states in a spin system 
              based on the number of spins and their respective spin quantum numbers. 
              The size of the state space is determined by the product of `(2 * spin + 1)` 
              for each spin in the system.
 @param total_spin_count The number of spins in the system.
 @param two_I An array containing `2 * I` values for each spin, where `I` is the spin quantum number.
 @return The total number of quantum states in the spin system as an integer.
 */
int number_of_states_(int total_spin_count, const int *two_I);

/*!
 @function get_single_spin_Ix_
 @abstract Creates the complex square matrix representation of the Ix operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Ix operator for the spin specified by `spin_index` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Ix operator will be stored.
 @param spin_index The index of the spin in the spin system for which the Ix operator is being calculated.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_index` is out of bounds, the function returns without performing any calculations.
 */

/*!
 @function get_single_spin_Ix_
 @abstract Creates the complex square matrix representation of the Ix operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Ix operator for the spin specified by `spin_index` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Ix operator will be stored.
 @param spin_index The index of the spin in the spin system for which the Ix operator is being calculated.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_index` is out of bounds, the function returns without performing any calculations.
 */
void get_single_spin_Ix_(double complex *operator, int spin_index, int *two_I, int total_spin_count);


/*!
 @function get_Ix_
@abstract Creates the complex square matrix representation of the total Ix operator for a subset of spins in a spin system.
@discussion This function generates the matrix representation of the Ix operator for a subset of spins specified by the `spin_indexes` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Ix operator will be stored.
 @param spin_indexes A integer array of the spin indexes in the spin system for which the Ix operator is being calculated.
 @param spin_count The total number of spins in the subsystem.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_indexes` are out of bounds, the function returns without performing any calculations.
 */
void get_Iz_(double complex *operator, int *spin_indexes, int spin_count, int *two_I, int total_spin_count);


/*!
 @function get_single_spin_Iy_
 @abstract Creates the complex square matrix representation of the Iy operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Iy operator for the spin specified by `spin_index` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Iy operator will be stored.
 @param spin_index The index of the spin in the spin system for which the Iy operator is being calculated.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_index` is out of bounds, the function returns without performing any calculations.
 */
void get_single_spin_Iy_(double complex *operator, int spin_index, int *two_I, int total_spin_count);

/*!
 @function get_Iy_
@abstract Creates the complex square matrix representation of the total Iy operator for a subset of spins in a spin system.
@discussion This function generates the matrix representation of the Iy operator for a subset of spins specified by the `spin_indexes` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Iy operator will be stored.
 @param spin_indexes A integer array of the spin indexes in the spin system for which the Iy operator is being calculated.
 @param spin_count The total number of spins in the subsystem.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_indexes` are out of bounds, the function returns without performing any calculations.
 */
void get_Iy_(double complex *operator, int *spin_indexes, int spin_count, int *two_I, int total_spin_count);



/*!
 @function get_single_spin_Iz_
 @abstract Creates the complex square matrix representation of the Iz operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Iz operator for the spin specified by `spin_index` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Iz operator will be stored.
 @param spin_index The index of the spin in the spin system for which the Iz operator is being calculated.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_index` is out of bounds, the function returns without performing any calculations.
 */
void get_single_spin_Iz_(double complex *operator, int spin_index, int *two_I, int total_spin_count);

/*!
 @function get_single_spin_Ip_
 @abstract Creates the complex square matrix representation of the Ip (I+) operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Ip operator for the spin specified by `spin_index` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Ip operator will be stored.
 @param spin_index The index of the spin in the spin system for which the Ip operator is being calculated.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_index` is out of bounds, the function returns without performing any calculations.
 */

/*!
 @function get_Iz_
@abstract Creates the complex square matrix representation of the total Iz operator for a subset of spins in a spin system.
@discussion This function generates the matrix representation of the Iz operator for a subset of spins specified by the `spin_indexes` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Iz operator will be stored.
 @param spin_indexes A integer array of the spin indexes in the spin system for which the Iz operator is being calculated.
 @param spin_count The total number of spins in the subsystem.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_indexes` are out of bounds, the function returns without performing any calculations.
 */
void get_Iz_(double complex *operator, int *spin_indexes, int spin_count, int *two_I, int total_spin_count);


void get_single_spin_Ip_(double complex *operator, int spin_index, int *two_I, int total_spin_count);

/*!
 @function get_single_spin_Im_
 @abstract Creates the complex square matrix representation of the Im (I−) operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Im operator for the spin specified by `spin_index` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Im operator will be stored.
 @param spin_index The index of the spin in the spin system for which the Im operator is being calculated.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_index` is out of bounds, the function returns without performing any calculations.
 */

/*!
 @function get_Ip_
@abstract Creates the complex square matrix representation of the total Ip operator for a subset of spins in a spin system.
@discussion This function generates the matrix representation of the Ip operator for a subset of spins specified by the `spin_indexes` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Ip operator will be stored.
 @param spin_indexes A integer array of the spin indexes in the spin system for which the Ip operator is being calculated.
 @param spin_count The total number of spins in the subsystem.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_indexes` are out of bounds, the function returns without performing any calculations.
 */
void get_Ip_(double complex *operator, int *spin_indexes, int spin_count, int *two_I, int total_spin_count);

void get_single_spin_Im_(double complex *operator, int spin_index, int *two_I, int total_spin_count);

/*!
 @function get_single_spin_Tlm_
 @abstract Creates the complex square matrix representation of the Tlm operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Tlm operator for the spin specified by `spin_index` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Tlm operator will be stored.
 @param spin_index The index of the spin in the spin system for which the Tlm operator is being calculated.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @param L The rank of the tensor operator.
 @param M The magnetic quantum number of the tensor operator.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_index` is out of bounds, the function returns without performing any calculations.
 */

/*!
 @function get_Im_
@abstract Creates the complex square matrix representation of the total Im operator for a subset of spins in a spin system.
@discussion This function generates the matrix representation of the Im operator for a subset of spins specified by the `spin_indexes` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Im operator will be stored.
 @param spin_indexes A integer array of the spin indexes in the spin system for which the Im operator is being calculated.
 @param spin_count The total number of spins in the subsystem.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_indexes` are out of bounds, the function returns without performing any calculations.
 */
void get_Ip_(double complex *operator, int *spin_indexes, int spin_count, int *two_I, int total_spin_count);


void get_single_spin_Tlm_(double complex *operator, int spin_index, int *two_I, int total_spin_count, int L, int M);
/*!
 @function get_single_spin_Tlm_unit_
 @abstract Creates the complex square matrix representation of the unit Tlm operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the unit Tlm operator for the spin specified by `spin_index` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients, normalization factors, and Kronecker delta products 
              to ensure proper coupling between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the unit Tlm operator will be stored.
 @param spin_index The index of the spin in the spin system for which the unit Tlm operator is being calculated.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @param L The rank of the tensor operator.
 @param M The magnetic quantum number of the tensor operator.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spin_index` is out of bounds, the function returns without performing any calculations.
 */
void get_single_spin_Tlm_unit_(double complex *operator, int spin_index, int *two_I, int total_spin_count, int L, int M);



/*!
 @function get_single_spin_C0_
 */
void get_single_spin_C0_(double complex *operator, int spin_index, int *two_I, int total_spin_count);

/*!
 @function get_single_spin_C2_
 */
void get_single_spin_C2_(double complex *operator, int spin_index, int *two_I, int total_spin_count);

/*!
 @function get_single_spin_C4_
 */
void get_single_spin_C4_(double complex *operator, int spin_index, int *two_I, int total_spin_count);


/*!
 @function get_Ef_
 @abstract Creates the complex square matrix representation of the identity operator for a fictitious spin-1/2 system.
 @discussion This function generates the matrix representation of the identity operator for a fictitious spin-1/2 system. 
              The operator acts on the specified states `r` and `s` in the spin system. The resulting matrix is stored 
              in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the identity operator will be stored.
 @param r The index of the first state.
 @param s The index of the second state.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 */
void get_Ef_(double complex *operator, int r, int s, int *two_I, int total_spin_count);

/*!
 @function get_Ixf_
 @abstract Creates the complex square matrix representation of the Ix operator for a fictitious spin-1/2 system.
 @discussion This function generates the matrix representation of the Ix operator for a fictitious spin-1/2 system. 
              The operator acts on the specified states `r` and `s` in the spin system. The resulting matrix is stored 
              in the provided `operator` array. The matrix elements are set to 0.5 for the off-diagonal elements 
              corresponding to the specified states and 0 for all other elements.
 @param operator A pointer to the array where the resulting complex square matrix for the Ix operator will be stored.
 @param r The index of the first state.
 @param s The index of the second state.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 */
void get_Ixf_(double complex *operator, int r, int s, int *two_I, int total_spin_count);

/*!
 @function get_Iyf_
 @abstract Creates the complex square matrix representation of the Iy operator for a fictitious spin-1/2 system.
 @discussion This function generates the matrix representation of the Iy operator for a fictitious spin-1/2 system. 
              The operator acts on the specified states `r` and `s` in the spin system. The resulting matrix is stored 
              in the provided `operator` array. The matrix elements are set to `0.5 * I` for the off-diagonal element 
              corresponding to `(r, s)`, `-0.5 * I` for `(s, r)`, and 0 for all other elements.
 @param operator A pointer to the array where the resulting complex square matrix for the Iy operator will be stored.
 @param r The index of the first state.
 @param s The index of the second state.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 */
void get_Iyf_(double complex *operator, int r, int s, int *two_I, int total_spin_count);

/*!
 @function get_Izf_
 @abstract Creates the complex square matrix representation of the Iz operator for a fictitious spin-1/2 system.
 @discussion This function generates the matrix representation of the Iz operator for a fictitious spin-1/2 system. 
              The operator acts on the specified states `r` and `s` in the spin system. The resulting matrix is stored 
              in the provided `operator` array. The matrix elements are set to `0.5` for the diagonal element 
              corresponding to state `s`, `-0.5` for state `r`, and 0 for all other elements.
 @param operator A pointer to the array where the resulting complex square matrix for the Iz operator will be stored.
 @param r The index of the first state.
 @param s The index of the second state.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 */
void get_Izf_(double complex *operator, int r, int s, int *two_I, int total_spin_count);

/*!
 @function get_Ipf_
 @abstract Creates the complex square matrix representation of the I+ (Iplus) operator for a fictitious spin-1/2 system.
 @discussion This function generates the matrix representation of the I+ operator for a fictitious spin-1/2 system. 
              The operator acts on the specified states `r` and `s` in the spin system. The resulting matrix is stored 
              in the provided `operator` array. The matrix element corresponding to `(s, r)` is set to 1, and all other 
              elements are set to 0.
 @param operator A pointer to the array where the resulting complex square matrix for the I+ operator will be stored.
 @param r The index of the first state.
 @param s The index of the second state.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 */
void get_Ipf_(double complex *operator, int r, int s, int *two_I, int total_spin_count);

/*!
 @function get_Imf_
 @abstract Creates the complex square matrix representation of the I− (Iminus) operator for a fictitious spin-1/2 system.
 @discussion This function generates the matrix representation of the I− operator for a fictitious spin-1/2 system. 
              The operator acts on the specified states `r` and `s` in the spin system. The resulting matrix is stored 
              in the provided `operator` array. The matrix element corresponding to `(r, s)` is set to 1, and all other 
              elements are set to 0.
 @param operator A pointer to the array where the resulting complex square matrix for the I− operator will be stored.
 @param r The index of the first state.
 @param s The index of the second state.
 @param two_I An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param total_spin_count The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 */
void get_Imf_(double complex *operator, int r, int s, int *two_I, int total_spin_count);

#endif // ANGULAR_MOMENTUM_H
