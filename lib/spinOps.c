#include "spinOps.h"

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

double Clebsch_(double j1,double m1,double j2,double m2,double j,double m)
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

double TLM_(double l,double m,double j1,double m1,double j2,double m2)
{
    double j;
    double element=0;
    if(j1==j2) {
        j = j1;
        double clebsch = Clebsch_(j,m2,l,m,j,m1);
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

double unitTLM_(double l,double m,double j1,double m1,double j2,double m2)
{
    double j;
    
    double element=0;
    if(j1==j2) {
        j = j1;
        element = Clebsch_(j2,m2,l,m,j1,m1)*sqrt(2*l+1)/sqrt(2*j+1);
    }
    return(element);
}

int numberOfStates_(int spinCount, int *spinsTimesTwo)
{
    /* Calculate size of state space */
    int nstates=1;
    for(unsigned int index = 0; index<spinCount; index++) {
        float spin = (float) spinsTimesTwo[index]/2.;
        nstates *= (unsigned int) (2. * spin + 1.);
    }
    return nstates;
}

float **createQuantumNumbers(int spinCount, int *spinsTimesTwo)
{
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);

    /* Create quantum numbers matrix */
    float **qnum = malloc(sizeof(float)*nstates*spinCount);
    
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
            }while(m <= spin);
        } while(state < nstates);
        x *= (2 * spin + 1.);
    }
    return qnum;
}

/* This routine calculates the product of delta(qnum[i][m1], qnum[i][m2]) for every spin i except spin iskip.*/

float systemDeltaProduct(float **qnum, int spinCount, int iskip, int bra, int ket)
{
    float delta=1.;
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
void getIx_(double complex **operator, int spinIndex, int *spinsTimesTwo, int spinCount)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    
    float **qnum = createQuantumNumbers(spinCount, spinsTimesTwo);
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);

    float spin = (float) spinsTimesTwo[spinIndex]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum, spinCount, spinIndex, bra, ket);
            if(del==0) operator[bra][ket] = 0;
            else {
                operator[bra][ket] = 1/ sqrt(2)*TLM_(1.,-1.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
                operator[bra][ket] -= 1/sqrt(2)*TLM_(1.,1.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
            }
        }
    }
    free(qnum);
}

/*!
 @function getIy
 @abstract create the Complex Square Matrix for Iy for the Spin in a Spin System
 @param spinIndex the index of spin in spin system.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Iy
 */
void getIy_(double complex **operator, int spinIndex, int *spinsTimesTwo, int spinCount)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    
    float **qnum = createQuantumNumbers(spinCount, spinsTimesTwo);
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);

    float spin = (float) spinsTimesTwo[spinIndex]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum, spinCount, spinIndex, bra, ket);
            if(del==0) operator[bra][ket] = 0;
            else {
                operator[bra][ket] = I/sqrt(2)*TLM_(1.,-1.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
                operator[bra][ket] += I/sqrt(2)*TLM_(1.,1.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
            }
        }
    }
    free(qnum);
}

/*!
 @function getIz
 @abstract create the Complex Square Matrix for Iz for the Spin in a Spin System
 @param spinIndex the index of spin in spin system.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Iz
 */
void getIz_(double complex **operator, int spinIndex, int *spinsTimesTwo, int spinCount)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    
    float **qnum = createQuantumNumbers(spinCount, spinsTimesTwo);
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);

    float spin = (float) spinsTimesTwo[spinIndex]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum, spinCount, spinIndex, bra, ket);
            if(del==0) operator[bra][ket] = 0;
            else operator[bra][ket] = TLM_(1.,0.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
        }
    }
    free(qnum);
}

/*!
 @function getIp
 @abstract create the Complex Square Matrix for Ip for the Spin in a Spin System
 @param spinIndex the index of spin in spin system.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Ip
 */
void getIp_(double complex **operator, int spinIndex, int *spinsTimesTwo, int spinCount)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    
    float **qnum = createQuantumNumbers(spinCount, spinsTimesTwo);
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);

    float spin = (float) spinsTimesTwo[spinIndex]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum, spinCount, spinIndex, bra, ket);
            if(del==0) operator[bra][ket] = 0;
            else operator[bra][ket] = - sqrt(2)*TLM_(1.,1.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
        }
    }
    free(qnum);
}

/*!
 @function getIm
 @abstract create the Complex Square Matrix for Im for the Spin in a Spin System
 @param spinIndex the index of spin in spin system.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Im
 */
void getIm_(double complex **operator, int spinIndex, int *spinsTimesTwo, int spinCount)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    
    float **qnum = createQuantumNumbers(spinCount, spinsTimesTwo);
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);

    float spin = (float) spinsTimesTwo[spinIndex]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum, spinCount, spinIndex, bra, ket);
            if(del==0) operator[bra][ket] = 0;
            else operator[bra][ket] = sqrt(2)*TLM_(1.,-1.,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
        }
    }
    free(qnum);
}


/*!
 @function getTlm
 @abstract create the Complex Square Matrix for Tlm for the Spin in a Spin System
 @param spinIndex the index of spin in spin system.
 @param spinsTimesTwo the integer array of 2*I for each spin in system.
 @param spinCount the count of spins in system.
 @result the Complex Square Matrix for Tlm
 */
void getTlm_(double complex **operator, int spinIndex, int *spinsTimesTwo, int spinCount, int L, int M)
{
    if(spinIndex<0 || spinIndex>spinCount-1) return; 
    
    float **qnum = createQuantumNumbers(spinCount, spinsTimesTwo);
    int nstates = numberOfStates_(spinCount, spinsTimesTwo);

    float spin = (float) spinsTimesTwo[spinIndex]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum, spinCount, spinIndex, bra, ket);
            if(del==0) operator[bra][ket] = 0;
            else operator[bra][ket] = TLM_(L,M,spin,qnum[spinIndex][bra],spin,qnum[spinIndex][ket]) * del;
        }
    }
    free(qnum);
}

