#include "spin.h"

/*!
 @function fac
 */
double fac(const double x) {
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
 */
double mypow(const double x, int n)
{
    double temp;
    if (n == 0) return (1.);
    temp = 1.;
    for (; n >= 1; n--) temp *= x;
    return (temp);
}

/*!
 @function deltaFunction
 */
float deltaFunction(const float m1, const float m2)
{
    float result = 1.;
    if (m1 != m2) result = 0.;
    return result;
}

/*!
 @function max
 */
double max(const double a, const double b, const double c)
{
    double m;
    if (a > b) m = a;
    else m = b;
    if (m < c) m = c;
    return (m);
}

/*!
 @function min
 */
double min(const double a, const double b, const double c)
{
    double m;
    if (a < b) m = a;
    else m = b;
    if (m > c) m = c;
    return (m);
}

/*!
 @function clebsch_
 */
double clebsch_(const double j1,const double m1,const double j2,const double m2,const double j,const double m)
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
 */
double tlm_(const double l,const double m,const double j1,const double m1,const double j2,const double m2)
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
 */
double unit_tlm_(const double l,const double m,const double j1,const double m1,const double j2,const double m2)
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
 @function number_of_states_
 */
int number_of_states_(int spin_count, int *i_times_2)
{
    /* Calculate size of state space */
    int nstates=1;
    for(int index = 0; index<spin_count; index++) {
        float spin = (float) i_times_2[index]/2.;
        nstates *= (unsigned int) (2. * spin + 1.);
    }
    return nstates;
}

/*!
 @function createQuantumNumbers
 */
float *createQuantumNumbers(int spin_count, int *i_times_2)
{
    int nstates = number_of_states_(spin_count, i_times_2);

    /* Create quantum numbers matrix */
    float *qnum_data = malloc(sizeof(float)*nstates*spin_count);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;

    double x = 1.;
    for(int index=0; index<spin_count; index++) {
        int state=0;
        float spin = (float) i_times_2[index]/2.;
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
 */
float systemDeltaProduct(float *qnum_data, int spin_count, int nstates, int iskip, int bra, int ket)
{
    float delta=1.;
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    for(int iSpin=0; iSpin<spin_count; iSpin++)
        if(iSpin!=iskip) delta *= deltaFunction(qnum[iSpin][bra], qnum[iSpin][ket]);
    return delta;
}

/*!
 @function get_single_spin_Ix_
 */
void get_single_spin_Ix_(double complex *operator, int spin_index, int *i_times_2, int spin_count)
{
    if(spin_index<0 || spin_index>spin_count-1) return; 
    int nstates = number_of_states_(spin_count, i_times_2);
    float *qnum_data = createQuantumNumbers(spin_count, i_times_2);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) i_times_2[spin_index]/2.;

    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spin_count, nstates, spin_index, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else {
                matrix[bra][ket] = 1/ sqrt(2)*tlm_(1.,-1.,spin,qnum[spin_index][bra],spin,qnum[spin_index][ket]) * del;
                matrix[bra][ket] -= 1/sqrt(2)*tlm_(1.,1.,spin,qnum[spin_index][bra],spin,qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}


/*!
 @function get_single_spin_Iy_
 */
void get_single_spin_Iy_(double complex *operator, int spin_index, int *i_times_2, int spin_count)
{
    if(spin_index<0 || spin_index>spin_count-1) return; 
    int nstates = number_of_states_(spin_count, i_times_2);
    float *qnum_data = createQuantumNumbers(spin_count, i_times_2);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) i_times_2[spin_index]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spin_count, nstates, spin_index, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else {
                matrix[bra][ket] = I/sqrt(2)*tlm_(1.,-1.,spin,qnum[spin_index][bra],spin,qnum[spin_index][ket]) * del;
                matrix[bra][ket] += I/sqrt(2)*tlm_(1.,1.,spin,qnum[spin_index][bra],spin,qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

/*!
 @function get_single_spin_Iz_
 */
void get_single_spin_Iz_(double complex *operator, int spin_index, int *i_times_2, int spin_count)
{
    if(spin_index<0 || spin_index>spin_count-1) return; 
    int nstates = number_of_states_(spin_count, i_times_2);
    float *qnum_data = createQuantumNumbers(spin_count, i_times_2);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) i_times_2[spin_index]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spin_count, nstates, spin_index, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else matrix[bra][ket] = tlm_(1.,0.,spin,qnum[spin_index][bra],spin,qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function get_single_spin_Ip_
 */
void get_single_spin_Ip_(double complex *operator, int spin_index, int *i_times_2, int spin_count)
{
    if(spin_index<0 || spin_index>spin_count-1) return; 
    int nstates = number_of_states_(spin_count, i_times_2);
    float *qnum_data = createQuantumNumbers(spin_count, i_times_2);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) i_times_2[spin_index]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spin_count, nstates, spin_index, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else matrix[bra][ket] = - sqrt(2)*tlm_(1.,1.,spin,qnum[spin_index][bra],spin,qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function get_single_spin_Im_
 */
void get_single_spin_Im_(double complex *operator, int spin_index, int *i_times_2, int spin_count)
{
    if(spin_index<0 || spin_index>spin_count-1) return; 
    int nstates = number_of_states_(spin_count, i_times_2);
    float *qnum_data = createQuantumNumbers(spin_count, i_times_2);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) i_times_2[spin_index]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spin_count, nstates, spin_index, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else matrix[bra][ket] = sqrt(2)*tlm_(1.,-1.,spin,qnum[spin_index][bra],spin,qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}


/*!
 @function get_single_spin_Tlm_
 */
void get_single_spin_Tlm_(double complex *operator, int spin_index, int *i_times_2, int spin_count, int L, int M)
{
    if(spin_index<0 || spin_index>spin_count-1) return; 
    int nstates = number_of_states_(spin_count, i_times_2);
    float *qnum_data = createQuantumNumbers(spin_count, i_times_2);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) i_times_2[spin_index]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spin_count, nstates, spin_index, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else matrix[bra][ket] = tlm_(L,M,spin,qnum[spin_index][bra],spin,qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function get_single_spin_Tlm_unit_
 */
void get_single_spin_Tlm_unit_(double complex *operator, int spin_index, int *i_times_2, int spin_count, int L, int M)
{
    if(spin_index<0 || spin_index>spin_count-1) return; 
    int nstates = number_of_states_(spin_count, i_times_2);
    float *qnum_data = createQuantumNumbers(spin_count, i_times_2);
    float (*qnum)[nstates] = (float (*)[nstates]) qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    float spin = (float) i_times_2[spin_index]/2.;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            float del = systemDeltaProduct(qnum_data, spin_count, nstates, spin_index, bra, ket);
            if(del==0) matrix[bra][ket] = 0;
            else matrix[bra][ket] = unit_tlm_(L,M,spin,qnum[spin_index][bra],spin,qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function getEf_
 */
void getEf_(double complex *operator, int r, int s, int *i_times_2, int spin_count)
{
    int nstates = number_of_states_(spin_count, i_times_2);
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
 */
void getIxf_(double complex *operator, int r, int s, int *i_times_2, int spin_count)
{
    int nstates = number_of_states_(spin_count, i_times_2);
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
 */
void getIyf_(double complex *operator, int r, int s, int *i_times_2, int spin_count)
{
    int nstates = number_of_states_(spin_count, i_times_2);
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
 */
void getIzf_(double complex *operator, int r, int s, int *i_times_2, int spin_count)
{
    int nstates = number_of_states_(spin_count, i_times_2);
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
 */
void getIpf_(double complex *operator, int r, int s, int *i_times_2, int spin_count)
{
    int nstates = number_of_states_(spin_count, i_times_2);
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
 */
void getImf_(double complex *operator, int r, int s, int *i_times_2, int spin_count)
{
    int nstates = number_of_states_(spin_count, i_times_2);
    double complex (*matrix)[nstates] = (double complex (*)[nstates]) operator;
    
    for(int bra=0; bra<nstates; bra++) {
        for(int ket=0; ket<nstates; ket++) {
            matrix[bra][ket] = 0;
			if((bra==r)&&(ket==s)) matrix[bra][ket] = 1;
        }
    }
}
