#ifndef SPIN_H
#define SPIN_H

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

/*!
 @function clebsch_
 @abstract Calculates the Clebsch-Gordon coefficients.
 @discussion This function computes the Clebsch-Gordon coefficients 
              `< j, m | j1, j2, m1, m2 >` using a routine adapted from the 
              Mathematica textbook (page 519). The Clebsch-Gordon coefficients 
              are used in quantum mechanics to describe the coupling of angular 
              momenta. The function ensures that the input values satisfy the 
              necessary conditions for valid coefficients.
 @param j1 The first angular momentum quantum number.
 @param m1 The magnetic quantum number associated with `j1`.
 @param j2 The second angular momentum quantum number.
 @param m2 The magnetic quantum number associated with `j2`.
 @param j The total angular momentum quantum number.
 @param m The total magnetic quantum number.
 @return The Clebsch-Gordon coefficient `< j, m | j1, j2, m1, m2 >` as a double. 
         Returns 0 if the input values do not satisfy the necessary conditions.
 */
double clebsch_(const double j1,const double m1,const double j2,const double m2,const double j,const double m);

/*!
 @function tlm_
 @abstract Evaluates the matrix element `<j1 m1|T(lm)|j2 m2>`.
 @discussion This function calculates the matrix element `<j1 m1|T(lm)|j2 m2>` 
              using the definition from Bowden and Hutchinson, J. Magn. Reson. 67, 403, 1986. 
              The calculation involves Clebsch-Gordon coefficients and reduced matrix elements. 
              The function assumes that `j1` equals `j2` for the calculation.
 @param l The rank of the tensor operator.
 @param m The order of the tensor operator.
 @param j1 The first angular momentum quantum number.
 @param m1 The magnetic quantum number associated with `j1`.
 @param j2 The second angular momentum quantum number.
 @param m2 The magnetic quantum number associated with `j2`.
 @return The matrix element `<j1 m1|T(lm)|j2 m2>` as a double. Returns 0 if `j1` is not equal to `j2`.
 */
double tlm_(double l,double m,double j1,double m1,double j2,double m2);

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
 @param j1 The first angular momentum quantum number.
 @param m1 The magnetic quantum number associated with `j1`.
 @param j2 The second angular momentum quantum number.
 @param m2 The magnetic quantum number associated with `j2`.
 @return The matrix element `<j1 m1|T_hat(lm)|j2 m2>` as a double. Returns 0 if `j1` is not equal to `j2`.
 */
double unit_tlm_(const double l,const double m,const double j1,const double m1,const double j2,const double m2);

/*!
 @function numberOfStates_
 @abstract Calculates the size of the state space for a spin system.
 @discussion This function computes the total number of quantum states in a spin system 
              based on the number of spins and their respective spin quantum numbers. 
              The size of the state space is determined by the product of `(2 * spin + 1)` 
              for each spin in the system.
 @param spinCount The number of spins in the system.
 @param spinsTimesTwo An array containing `2 * I` values for each spin, where `I` is the spin quantum number.
 @return The total number of quantum states in the spin system as an integer.
 */
int numberOfStates_(int spinCount, int *spinsTimesTwo);

/*!
 @function getIx_
 @abstract Creates the complex square matrix representation of the Ix operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Ix operator for the spin specified by `spinIndex` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Ix operator will be stored.
 @param spinIndex The index of the spin in the spin system for which the Ix operator is being calculated.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spinIndex` is out of bounds, the function returns without performing any calculations.
 */

/*!
 @function getIx_
 @abstract Creates the complex square matrix representation of the Ix operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Ix operator for the spin specified by `spinIndex` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Ix operator will be stored.
 @param spinIndex The index of the spin in the spin system for which the Ix operator is being calculated.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spinIndex` is out of bounds, the function returns without performing any calculations.
 */
void getIx_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount);

/*!
 @function getIy_
 @abstract Creates the complex square matrix representation of the Iy operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Iy operator for the spin specified by `spinIndex` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Iy operator will be stored.
 @param spinIndex The index of the spin in the spin system for which the Iy operator is being calculated.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spinIndex` is out of bounds, the function returns without performing any calculations.
 */
void getIy_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount);

/*!
 @function getIz_
 @abstract Creates the complex square matrix representation of the Iz operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Iz operator for the spin specified by `spinIndex` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Iz operator will be stored.
 @param spinIndex The index of the spin in the spin system for which the Iz operator is being calculated.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spinIndex` is out of bounds, the function returns without performing any calculations.
 */
void getIz_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount);

/*!
 @function getIp_
 @abstract Creates the complex square matrix representation of the Ip (I+) operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Ip operator for the spin specified by `spinIndex` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Ip operator will be stored.
 @param spinIndex The index of the spin in the spin system for which the Ip operator is being calculated.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spinIndex` is out of bounds, the function returns without performing any calculations.
 */
void getIp_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount);

/*!
 @function getIm_
 @abstract Creates the complex square matrix representation of the Im (I−) operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Im operator for the spin specified by `spinIndex` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Im operator will be stored.
 @param spinIndex The index of the spin in the spin system for which the Im operator is being calculated.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spinIndex` is out of bounds, the function returns without performing any calculations.
 */
void getIm_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount);

/*!
 @function getTlm_
 @abstract Creates the complex square matrix representation of the Tlm operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the Tlm operator for the spin specified by `spinIndex` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients and Kronecker delta products to ensure proper coupling 
              between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the Tlm operator will be stored.
 @param spinIndex The index of the spin in the spin system for which the Tlm operator is being calculated.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @param L The rank of the tensor operator.
 @param M The magnetic quantum number of the tensor operator.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spinIndex` is out of bounds, the function returns without performing any calculations.
 */
void getTlm_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount, int L, int M);
/*!
 @function getTlm_unit_
 @abstract Creates the complex square matrix representation of the unit Tlm operator for a single spin in a spin system.
 @discussion This function generates the matrix representation of the unit Tlm operator for the spin specified by `spinIndex` 
              in a spin system. The matrix is constructed in the basis of quantum states for the system, and the 
              calculation involves Clebsch-Gordon coefficients, normalization factors, and Kronecker delta products 
              to ensure proper coupling between states. The resulting matrix is stored in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the unit Tlm operator will be stored.
 @param spinIndex The index of the spin in the spin system for which the unit Tlm operator is being calculated.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @param L The rank of the tensor operator.
 @param M The magnetic quantum number of the tensor operator.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 @note If `spinIndex` is out of bounds, the function returns without performing any calculations.
 */
void getTlm_unit_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount, int L, int M);

/*!
 @function getEf_
 @abstract Creates the complex square matrix representation of the identity operator for a fictitious spin-1/2 system.
 @discussion This function generates the matrix representation of the identity operator for a fictitious spin-1/2 system. 
              The operator acts on the specified states `r` and `s` in the spin system. The resulting matrix is stored 
              in the provided `operator` array.
 @param operator A pointer to the array where the resulting complex square matrix for the identity operator will be stored.
 @param r The index of the first state.
 @param s The index of the second state.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 */
void getEf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);

/*!
 @function getIxf_
 @abstract Creates the complex square matrix representation of the Ix operator for a fictitious spin-1/2 system.
 @discussion This function generates the matrix representation of the Ix operator for a fictitious spin-1/2 system. 
              The operator acts on the specified states `r` and `s` in the spin system. The resulting matrix is stored 
              in the provided `operator` array. The matrix elements are set to 0.5 for the off-diagonal elements 
              corresponding to the specified states and 0 for all other elements.
 @param operator A pointer to the array where the resulting complex square matrix for the Ix operator will be stored.
 @param r The index of the first state.
 @param s The index of the second state.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 */
void getIxf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);

/*!
 @function getIyf_
 @abstract Creates the complex square matrix representation of the Iy operator for a fictitious spin-1/2 system.
 @discussion This function generates the matrix representation of the Iy operator for a fictitious spin-1/2 system. 
              The operator acts on the specified states `r` and `s` in the spin system. The resulting matrix is stored 
              in the provided `operator` array. The matrix elements are set to `0.5 * I` for the off-diagonal element 
              corresponding to `(r, s)`, `-0.5 * I` for `(s, r)`, and 0 for all other elements.
 @param operator A pointer to the array where the resulting complex square matrix for the Iy operator will be stored.
 @param r The index of the first state.
 @param s The index of the second state.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 */
void getIyf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);

/*!
 @function getIzf_
 @abstract Creates the complex square matrix representation of the Iz operator for a fictitious spin-1/2 system.
 @discussion This function generates the matrix representation of the Iz operator for a fictitious spin-1/2 system. 
              The operator acts on the specified states `r` and `s` in the spin system. The resulting matrix is stored 
              in the provided `operator` array. The matrix elements are set to `0.5` for the diagonal element 
              corresponding to state `s`, `-0.5` for state `r`, and 0 for all other elements.
 @param operator A pointer to the array where the resulting complex square matrix for the Iz operator will be stored.
 @param r The index of the first state.
 @param s The index of the second state.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 */
void getIzf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);

/*!
 @function getIpf_
 @abstract Creates the complex square matrix representation of the I+ (Iplus) operator for a fictitious spin-1/2 system.
 @discussion This function generates the matrix representation of the I+ operator for a fictitious spin-1/2 system. 
              The operator acts on the specified states `r` and `s` in the spin system. The resulting matrix is stored 
              in the provided `operator` array. The matrix element corresponding to `(s, r)` is set to 1, and all other 
              elements are set to 0.
 @param operator A pointer to the array where the resulting complex square matrix for the I+ operator will be stored.
 @param r The index of the first state.
 @param s The index of the second state.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 */
void getIpf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);

/*!
 @function getImf_
 @abstract Creates the complex square matrix representation of the I− (Iminus) operator for a fictitious spin-1/2 system.
 @discussion This function generates the matrix representation of the I− operator for a fictitious spin-1/2 system. 
              The operator acts on the specified states `r` and `s` in the spin system. The resulting matrix is stored 
              in the provided `operator` array. The matrix element corresponding to `(r, s)` is set to 1, and all other 
              elements are set to 0.
 @param operator A pointer to the array where the resulting complex square matrix for the I− operator will be stored.
 @param r The index of the first state.
 @param s The index of the second state.
 @param spinsTimesTwo An array containing `2 * I` values for each spin in the system, where `I` is the spin quantum number.
 @param spinCount The total number of spins in the system.
 @return This function does not return a value. The resulting matrix is stored in the `operator` array.
 */
void getImf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);


/*!
 @function mypow
 @abstract Computes the power of a number.
 @discussion This function calculates `x` raised to the power of `n`. 
              If `n` is 0, the function returns 1. For positive values of `n`, 
              the function iteratively multiplies `x` by itself `n` times.
 @param x The base value.
 @param n The exponent (non-negative integer).
 @return The result of `x` raised to the power of `n` as a double.
 */
double mypow(const double x, int n);

/*!
 @function fac
 @abstract Computes the factorial of a non-negative number.
 @discussion This function calculates the factorial of a non-negative number `x`. 
              If `x` is not an integer, it is truncated to its integer part. 
              If `x` is negative, an error message is printed, and the function 
              returns 0. For `x = 0`, the function returns 1 (by definition).
 @param x The input number for which the factorial is to be computed.
 @return The factorial of the input number `x` as a double. Returns 0 if `x` is negative.
 */
double fac(const double x);

#ifdef __cplusplus
}
#endif

#endif // SPIN_H
