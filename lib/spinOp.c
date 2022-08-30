#include <stdio.h>
#include <math.h>

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

double Clebsch(double j1,double m1,double j2,double m2,double j,double m)
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

double TLM(double l,double m,double j1,double m1,double j2,double m2)
{
    double j;
    double element=0;
    if(j1==j2) {
        j = j1;
        double clebsch = Clebsch(j,m2,l,m,j,m1);
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

double unitTLM(double l,double m,double j1,double m1,double j2,double m2)
{
    double j;
    
    double element=0;
    if(j1==j2) {
        j = j1;
        element = Clebsch(j2,m2,l,m,j1,m1)*sqrt(2*l+1)/sqrt(2*j+1);
    }
    return(element);
}

