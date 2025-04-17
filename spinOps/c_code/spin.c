#include "spin.h"

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
double fac(double x) {
    // Validate input
    if (x < 0) {
        fprintf(stderr, "Error: illegal argument x = %g in factorial. Factorial is undefined for negative numbers.\n", x);
        return 0;
    }

    // Handle edge case for x = 0
    if (x == 0) {
        return 1.0;
    }

    // Truncate x to its integer part
    int ix = (int)x;

    // Compute factorial iteratively
    double result = 1.0;
    for (int i = 2; i <= ix; i++) {
        result *= i;
    }

    return result;
}

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
double mypow(double x, int n)
{
    double temp;
    if (n == 0) return (1.);
    temp = 1.;
    for (; n >= 1; n--) temp *= x;
    return (temp);
}

/*!
 @function deltaFunction
 @abstract Computes the Kronecker delta for two floating-point numbers.
 @discussion This function evaluates the Kronecker delta, which is 1 if the two 
              input values `m1` and `m2` are equal, and 0 otherwise.
 @param m1 The first input value.
 @param m2 The second input value.
 @return Returns 1 if `m1` equals `m2`, otherwise returns 0.
 */
float deltaFunction(float m1, float m2)
{
    float result = 1.;
    if (m1 != m2) result = 0.;
    return result;
}

/*!
 @function max
 @abstract Computes the maximum of three double values.
 @discussion This function compares three double values `a`, `b`, and `c` 
              and returns the largest of the three.
 @param a The first value.
 @param b The second value.
 @param c The third value.
 @return The largest value among `a`, `b`, and `c`.
 */
double max(double a, double b, double c)
{
    double m;
    if (a > b) m = a;
    else m = b;
    if (m < c) m = c;
    return (m);
}

/*!
 @function min
 @abstract Computes the minimum of three double values.
 @discussion This function compares three double values `a`, `b`, and `c` 
              and returns the smallest of the three.
 @param a The first value.
 @param b The second value.
 @param c The third value.
 @return The smallest value among `a`, `b`, and `c`.
 */
double min(double a, double b, double c)
{
    double m;
    if (a < b) m = a;
    else m = b;
    if (m > c) m = c;
    return (m);
}

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
double clebsch_(double j1,double m1,double j2,double m2,double j,double m)
{
    double C1 = 0.0, C2, C3, temp;
    double cg = 0.0;
    int imin, imax, k;
    
    if(fabs(m) > j) return(0.);
    if(m1+m2 == m) {
        imin = (int) max(0., j2-j-m1, j1-j+m2);
        imax = (int) min(j1+j2-j, j1-m1, j2+m2);
        for(k=imin; k<=imax; k++) {
            temp = fac((double)k) * fac(j1 + j2 - j - (double)k)
            * fac(j1 - m1 - (double)k) * fac( j2 + m2 - (double)k)
            * fac(j - j2 + m1 + (double)k) * fac(j - j1 - m2 + (double)k);
            C1 += pow(-1, k) / temp;
        }
        C2 = fac(-j+j1+j2) * fac(j-j1+j2) * fac(j+j1-j2) * (2*j+1) / fac(1.+j+j1+j2);
        C3 = fac(j-m) * fac(j+m) * fac(j1-m1) * fac(j1+m1) * fac(j2-m2) * fac(j2+m2);
        cg = C1 * sqrt(C2 * C3);
    }
    return(cg);
}


/*!
 @function tlm_
 @abstract Evaluates the matrix element `<j1 m1|T(lm)|j2 m2>`.
 @discussion This function calculates the matrix element `<j1 m1|T(lm)|j2 m2>` 
              using the definition from Bowden and Hutchinson, J. Magn. Reson. 67, 403, 1986. 
              The calculation involves Clebsch-Gordon coefficients and reduced matrix elements. 
              The function assumes that `j1` equals `j2` for the calculation.
 @param l The rank of the tensor operator.
 @param m The magnetic quantum number of the tensor operator.
 @param j1 The first angular momentum quantum number.
 @param m1 The magnetic quantum number associated with `j1`.
 @param j2 The second angular momentum quantum number.
 @param m2 The magnetic quantum number associated with `j2`.
 @return The matrix element `<j1 m1|T(lm)|j2 m2>` as a double. Returns 0 if `j1` is not equal to `j2`.
 */
double tlm_(double l,double m,double j1,double m1,double j2,double m2)
{
    double j;
    double element=0;
    if(j1==j2) {
        j = j1;
        double clebsch = clebsch_(j,m2,l,m,j,m1);
        if(clebsch!=0.0) {
            double rme = fac(l) * fac(l) * fac(2*j+l+1);
            rme /= pow(2.,l) * fac(2*l) * fac(2*j - l);
            rme = sqrt(rme);
            element = clebsch * rme / sqrt(2*j+1);
        }
    }
    return(element);
}

/*!
 @function unit_tlm_
 @abstract Evaluates the matrix element `<j1 m1|T_hat(lm)|j2 m2>` for unit tensors.
 @discussion This function calculates the matrix element `<j1 m1|T_hat(lm)|j2 m2>` 
              using the definition of unit tensors from Bowden and Hutchinson, 
              J. Magn. Reson. 67, 403, 1986. The calculation involves Clebsch-Gordon 
              coefficients and normalization factors. The function assumes that 
              `j1` equals `j2` for the calculation.
 @param l The rank of the tensor operator.
 @param m The magnetic quantum number of the tensor operator.
 @param j1 The first angular momentum quantum number.
 @param m1 The magnetic quantum number associated with `j1`.
 @param j2 The second angular momentum quantum number.
 @param m2 The magnetic quantum number associated with `j2`.
 @return The matrix element `<j1 m1|T_hat(lm)|j2 m2>` as a double. Returns 0 if `j1` is not equal to `j2`.
 */
double unit_tlm_(double l,double m,double j1,double m1,double j2,double m2)
{
    double j;
    
    double element=0;
    if(j1==j2) {
        j = j1;
        element = clebsch_(j2,m2,l,m,j1,m1)*sqrt(2*l+1)/sqrt(2*j+1);
    }
    return(element);
}

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
int numberOfStates_(int spinCount, int *spinsTimesTwo)
{
    /* Calculate size of state space */
    int nstates=1;
    for(int index = 0; index<spinCount; index++) {
        float spin = (float) spinsTimesTwo[index]/2.;
        nstates *= (unsigned int) (2. * spin + 1.);
    }
    return nstates;
}

/*!
 @function createQuantumNumbers
 @abstract Creates the quantum numbers matrix for a spin system.
 @discussion This function generates a matrix of quantum numbers for a spin system, 
              where each row corresponds to a spin and each column corresponds to a 
              quantum state. The quantum numbers are calculated based on the spin 
              quantum numbers provided in the `spinsTimesTwo` array. The matrix is 
              stored in a dynamically allocated array.
 @param spinCount The number of spins in the system.
 @param spinsTimesTwo An array containing `2 * I` values for each spin, where `I` is the spin quantum number.
 @return A pointer to the dynamically allocated array containing the quantum numbers matrix.
         The caller is responsible for freeing the allocated memory.
 */
float *createQuantumNumbers(int spinCount, int *spinsTimesTwo)
{
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);

    /* Create quantum numbers matrix */
    float *qnum_data = malloc(sizeof(float)*nstates*spinCount);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;

    double x = 1.;
    for(int index=0; index<spinCount; index++) {
        int state=0;
        float spin = (float) spinsTimesTwo[index]/2.;
        do {float m = - spin;
            do {
                qnum[index][state] = (float) m;
                state++;
                double ip;
                if(modf( (double) state/x,&ip) == 0.) m++;
            } while(m <= spin);
        } while(state < nstates);
        x *= (2 * spin + 1.);
    }
    return qnum_data;
}

/*!
 @function systemDeltaProduct
 @abstract Calculates the product of Kronecker deltas for a spin system.
 @discussion This function computes the product of Kronecker delta values 
              `delta(qnum[i][bra], qnum[i][ket])` for all spins `i` in the system, 
              excluding the spin specified by `iskip`. The Kronecker delta is 1 
              if the two quantum numbers are equal and 0 otherwise.
 @param qnum_data A pointer to the quantum numbers matrix for the spin system.
 @param spinCount The total number of spins in the system.
 @param nstates The total number of quantum states in the system.
 @param iskip The index of the spin to exclude from the calculation.
 @param bra The index of the bra state.
 @param ket The index of the ket state.
 @return The product of Kronecker delta values as a float.
 */
float systemDeltaProduct(float *qnum_data, int spinCount, int nstates, int iskip, int bra, int ket)
{
    float delta=1.;
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    for(int iSpin=0; iSpin<spinCount; iSpin++)
        if(iSpin!=iskip) delta *= deltaFunction(qnum[iSpin][bra], qnum[iSpin][ket]);
    return delta;
}

/*!
 @function getIx_
 @abstract Creates the complex square matrix representation of the Ix operator for a specific spin in a spin system.
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
void getIx_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    float *qnum_data = createQuantumNumbers(spinCount, spinsTimesTwo);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) spinsTimesTwo[spinIndex]/2.;

    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spinCount, nstates, spinIndex, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else {
                matrix[bra][ket] = 1/ sqrt(2)*tlm_(1.,-1.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
                matrix[bra][ket] -= 1/sqrt(2)*tlm_(1.,1.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
            }
        }
    }
    free(qnum_data);
}


/*!
 @function getIy_
 @abstract Creates the complex square matrix representation of the Iy operator for a specific spin in a spin system.
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
void getIy_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    float *qnum_data = createQuantumNumbers(spinCount, spinsTimesTwo);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) spinsTimesTwo[spinIndex]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spinCount, nstates, spinIndex, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else {
                matrix[bra][ket] = I/sqrt(2)*tlm_(1.,-1.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
                matrix[bra][ket] += I/sqrt(2)*tlm_(1.,1.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

/*!
 @function getIz_
 @abstract Creates the complex square matrix representation of the Iz operator for a specific spin in a spin system.
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
void getIz_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    float *qnum_data = createQuantumNumbers(spinCount, spinsTimesTwo);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) spinsTimesTwo[spinIndex]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spinCount, nstates, spinIndex, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else matrix[bra][ket] = tlm_(1.,0.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function getIp_
 @abstract Creates the complex square matrix representation of the Ip (I+) operator for a specific spin in a spin system.
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
void getIp_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    float *qnum_data = createQuantumNumbers(spinCount, spinsTimesTwo);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) spinsTimesTwo[spinIndex]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spinCount, nstates, spinIndex, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else matrix[bra][ket] = - sqrt(2)*tlm_(1.,1.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function getIm_
 @abstract Creates the complex square matrix representation of the Im (I−) operator for a specific spin in a spin system.
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
void getIm_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    float *qnum_data = createQuantumNumbers(spinCount, spinsTimesTwo);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) spinsTimesTwo[spinIndex]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spinCount, nstates, spinIndex, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else matrix[bra][ket] = sqrt(2)*tlm_(1.,-1.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
        }
    }
    free(qnum_data);
}


/*!
 @function getTlm_
 @abstract Creates the complex square matrix representation of the Tlm operator for a specific spin in a spin system.
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
void getTlm_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount, int L, int M)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    float *qnum_data = createQuantumNumbers(spinCount, spinsTimesTwo);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) spinsTimesTwo[spinIndex]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spinCount, nstates, spinIndex, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else matrix[bra][ket] = tlm_(L,M,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function getTlm_unit_
 @abstract Creates the complex square matrix representation of the unit Tlm operator for a specific spin in a spin system.
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
void getTlm_unit_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount, int L, int M)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    float *qnum_data = createQuantumNumbers(spinCount, spinsTimesTwo);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) spinsTimesTwo[spinIndex]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spinCount, nstates, spinIndex, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else matrix[bra][ket] = unit_tlm_(L,M,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
        }
    }
    free(qnum_data);
}

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
void getEf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount)
{
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            matrix[bra][ket] = 0;
			if(bra==ket&&ket==s) matrix[bra][ket] = 1;
			else if(bra==ket&&ket==r) matrix[bra][ket] = 1;
        }
    }
}

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
void getIxf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount)
{
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            matrix[bra][ket] = 0;
			if((bra==r)&&(ket==s)) matrix[bra][ket] = .5;
			else if((bra==s)&&(ket==r)) matrix[bra][ket] = .5;
        }
    }
}

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
void getIyf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount)
{
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            matrix[bra][ket] = 0;
			if(bra==r&&ket==s) matrix[bra][ket] = .5*I;
			else if(bra==s&&ket==r) matrix[bra][ket] = -.5*I;
        }
    }
}

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
void getIzf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount)
{
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            matrix[bra][ket] = 0;
			if(bra==ket&&ket==s) matrix[bra][ket] = .5;
			else if(bra==ket&&ket==r) matrix[bra][ket] = -.5;
        }
    }
}

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
void getIpf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount)
{
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            matrix[bra][ket] = 0;
			if((ket==r)&&(bra==s)) matrix[bra][ket] = 1;
        }
    }
}

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
void getImf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount)
{
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            matrix[bra][ket] = 0;
			if((bra==r)&&(ket==s)) matrix[bra][ket] = 1;
        }
    }
}
