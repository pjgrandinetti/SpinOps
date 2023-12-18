#include "spin.h"

double fac(double x)
{
    if (x < 0) {
        fprintf(stderr, "illegal argument x = %g in factorial...\n",x);
        return 0;
    }
    int ix = (int) x;
    double sum = 1;
    for (; ix > 1; ix--) sum *= ix;
    return sum;
}

/* power function */

double mypow(double x, int n)
{
    double temp;
    if(n==0) return(1.);
    temp = 1.;
    for (; n >= 1; n--) temp *= x;
    return(temp);
}

float deltaFunction(float m1,float m2)
{
    float result=1.;
    if(m1!=m2) result=0.;
    return result;
}


double max(double a,double b,double c)
{
    double m;
    if(a > b) m = a;
    else m = b;
    if(m < c) m = c;
    return(m);
}

double min(double a,double b,double c)
{
    double m;
    if(a < b) m = a;
    else m = b;
    if(m > c) m = c;
    return(m);
}

/* This routine calculates the Clebsch-Gordon coefficients */
/* < j, m | j1, j2, m1, m2>, using a routine taken from the */
/* Mathematica textbook, page 519. */

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


/* This routines evaluates the matrix element <j1 m1|T(lm)|j2 m2> using */
/* definition from Bowden and Hutchinson, J. magn. reson. 67, 403, 1986 */

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

/* This routines evaluates the matrix element <j1 m1|T_hat(lm)|j2 m2> using */
/* definition of unit Tensors from Bowden and Hutchinson, J. magn. reson. 67, 403, 1986 */

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

/* This routine calculates the product of delta(qnum[i][m1], qnum[i][m2]) for every spin i except spin iskip.*/

float systemDeltaProduct(float *qnum_data, int spinCount, int nstates, int iskip, int bra, int ket)
{
    float delta=1.;
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    for(int iSpin=0; iSpin<spinCount; iSpin++)
        if(iSpin!=iskip) delta *= deltaFunction(qnum[iSpin][bra], qnum[iSpin][ket]);
    return delta;
}

/*!
 @function getIx
 @abstract create the Complex Square Matrix for Ix for the Spin in a Spin System
 @param spinIndex the index of spin in spin system.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Ix
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
 @function getIy
 @abstract create the Complex Square Matrix for Iy for the Spin in a Spin System
 @param spinIndex the index of spin in spin system.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Iy
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
 @function getIz
 @abstract create the Complex Square Matrix for Iz for the Spin in a Spin System
 @param spinIndex the index of spin in spin system.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Iz
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
 @function getIp
 @abstract create the Complex Square Matrix for Ip for the Spin in a Spin System
 @param spinIndex the index of spin in spin system.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Ip
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
 @function getIm
 @abstract create the Complex Square Matrix for Im for the Spin in a Spin System
 @param spinIndex the index of spin in spin system.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Im
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
 @function getTlm
 @abstract create the Complex Square Matrix for Tlm for the Spin in a Spin System
 @param spinIndex the index of spin in spin system.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Tlm
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
 @function getTlm_unit
 @abstract create the Complex Square Matrix for unit Tlm for the Spin in a Spin System
 @param spinIndex the index of spin in spin system.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Tlm
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
 @abstract create the Complex Square Matrix for Fictitious Spin 1/2 
 Identity Operator in a Spin System
 @param r state index.
 @param s state index.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Im
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
 @abstract create the Complex Square Matrix for Fictitious Spin 1/2 
 Ix Operator in a Spin System
 @param r state index.
 @param s state index.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Im
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
 @abstract create the Complex Square Matrix for Fictitious Spin 1/2 
 Iy Operator in a Spin System
 @param r state index.
 @param s state index.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Im
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
 @abstract create the Complex Square Matrix for Fictitious Spin 1/2 
 Iz Operator in a Spin System
 @param r state index.
 @param s state index.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Im
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
 @abstract create the Complex Square Matrix for Fictitious Spin 1/2 
 Iplus Operator in a Spin System
 @param r state index.
 @param s state index.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Im
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
 @abstract create the Complex Square Matrix for Fictitious Spin 1/2 
 Iminus Operator in a Spin System
 @param r state index.
 @param s state index.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Im
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
