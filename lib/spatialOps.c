#include "spatialOps.h"

/*!
 @function getrho2_pas
 @abstract create the Complex vector for the spherical tensor of rank 2
 @param zeta the traceless 2nd-rank symetric tensor anisotropy.
 @param eta the traceless 2nd-rank symetric tensor asymmetrt.
 */
void getrho1_pas_(double complex *tensor, double zeta, double eta)
{
	double complex *rho	= tensor + 1;
	for (int k = -1; k <=1; k++) rho[k] = 0.;
	rho[0] =  zeta;   
	rho[1] = zeta*eta;
}

/*!
 @function getrho2_pas
 @abstract create the Complex vector for the spherical tensor of rank 2
 @param zeta the traceless 2nd-rank symetric tensor anisotropy.
 @param eta the traceless 2nd-rank symetric tensor asymmetrt.
 */
void getrho2_pas_(double complex *tensor, double zeta, double eta)
{
	double complex *rho	= tensor + 2;
	for (int k = -2; k <= 2; k++) rho[k] = 0.;
	rho[0] = zeta/sqrt(6.);   
	rho[2] = rho[-2] = rho[0] * eta/sqrt(6.);
}

/* calculate Wigner rotation matrices */

double wigner_d_(double l,double m1,double m2,double beta)
{
	if(l==2) {
		if(m1==2) {
			if(m2==2) {
				double cx = cos(beta);
				return( (1+cx)*(1.+cx)/4.);
				}
			else if(m2==1) {
				double sx = sin(beta);
				double cx = cos(beta);
				return( -sx*(1.+cx)/2.);
				}
			else if(m2==0) {
				double sx = sin(beta);
				return( 0.6123724355*sx*sx);
				}
			else if(m2==-1) {
				double sx = sin(beta);
				double cx = cos(beta);
				return( -sx*(1.-cx)/2.);
				}
			else if(m2==-2) {
				double cx = cos(beta);
				return( (1-cx)*(1.-cx)/4.);
				}
			}
		else if(m1==-2) {
			if(m2==2) {
				double cx = cos(beta);
				return( (1-cx)*(1.-cx)/4.);
				}
			else if(m2==1) {
				double sx = sin(beta);
				double cx = cos(beta);
				return(sx*(1.-cx)/2.);
				}
			else if(m2==0) {
				double sx = sin(beta);
				return( 0.6123724355*sx*sx);
				}
			else if(m2==-1) {
				double sx = sin(beta);
				double cx = cos(beta);
				return(sx*(1.+cx)/2.);
				}
			else if(m2==-2) {
				double cx = cos(beta);
				return( (1+cx)*(1.+cx)/4.);
				}
			}
		else if(m1==1) {
			if(m2==2) {
				double sx = sin(beta);
				double cx = cos(beta);
				return( sx*(1+cx)/2.);
				}
			else if(m2==1) {
				double cx = cos(beta);
				return((2*cx*cx+cx-1.)/2.);
				}
			else if(m2==0) {
				double sx = sin(beta);
				double cx = cos(beta);
				return(-1.224744871*sx*cx);
				}
			else if(m2==-1) {
				double cx = cos(beta);
				return(-(2*cx*cx-cx-1.)/2.);
				}
			else if(m2==-2) {
				double sx = sin(beta);
				double cx = cos(beta);
				return( -sx*(1-cx)/2.);
				}
			}
		else if(m1==0) {
			if(m2==2) {
				double sx = sin(beta);
				return(0.6123724355*sx*sx);
				}
			else if(m2==1) {
				double sx = sin(beta);
				double cx = cos(beta);
				return(1.224744871*sx*cx);
				}
			else if(m2==0) {
				double cx = cos(beta);
				return(1.5*cx*cx- .5);
				}
			else if(m2==-1) {
				double sx = sin(beta);
				double cx = cos(beta);
				return(-1.224744871*sx*cx);
				}
			else if(m2==-2) {
				double sx = sin(beta);
				return(0.6123724355*sx*sx);
				}
			}
		else if(m1==-1) {
			if(m2==2) {
				double sx = sin(beta);
				double cx = cos(beta);
				return( sx*(1-cx)/2.);
				}
			else if(m2==1) {
				double cx = cos(beta);
				return(-(2*cx*cx-cx-1.)/2.);
				}
			else if(m2==0) {
				double sx = sin(beta);
				double cx = cos(beta);
				return(1.224744871*sx*cx);
				}
			else if(m2==-1) {
				double cx = cos(beta);
				return((2*cx*cx+cx-1.)/2.);
				}
			else if(m2==-2) {
				double sx = sin(beta);
				double cx = cos(beta);
				return( -sx*(1+cx)/2.);
				}
			}
		}
	else {
		double fac(double x), mypow(double x, int n);

		double sx = sin(beta/2.);
		double cx = cos(beta/2.);
		double sum = 0.;
		int sign = 1;
	                          
		for (int k = 0; k <= l - m1; k++) {
			double k1 = (int)(l - m1 - k);
			double k2 = (int)(l + m2 - k);
			double k3 = (int)(k + m1 - m2);

			if ( k1 >= 0 && k2 >= 0 && k3 >= 0) {   
			 	int n1 = (int)(2 * l + m2 - m1 - 2 * k);
				int n2 = (int)(m1 - m2 + 2 * k);
				double x = mypow(cx, n1);
				double y = mypow(sx, n2);
				sum += sign * x * y / (fac((double)k1) * fac((double)k2) * fac((double)k3) * fac((double)k)); 
				}
		sign = -sign;
		}
		double f = fac(l+m2) * fac(l-m2) * fac(l+m1) * fac(l-m1);
		f = sqrt(f);
		return(sum * f);
		}
	return(0);
}

double complex DLM_(double l,double  m1,double m2, double alpha, double beta, double gamma)
{
	double pha = m1 * alpha + m2 *gamma;
	double db = wigner_d_(l, m1, m2, beta);                                   
	return cos(pha) * db -I* sin(pha) * db;
}

/* Rotational transformation from one frame to another frame */

void Rot_(double j, double complex *initial, double alpha, double beta, double gamma, double complex *final)
{
	int m1, m2;
	double complex d;

	if((alpha==0.)&&(beta==0.)&&(gamma==0.)) {
		for (m2 = -(int)j; m2 <= (int)j; m2++) final[m2] = initial[m2];
		return;
		} 
	
	if(j==2) {
		double pha, db;
		for (m2 = -2; m2 <= 2; m2++) {
			final[m2] =0.;
			for (m1 = -2; m1 <= 2; m1++) {
				db = wigner_d_(j, m1, m2, beta);                                   
				pha = m1 * alpha + m2 * gamma;
				d = cos(pha) * db - I* sin(pha) * db;
	        	final[m2] += d * initial[m1];
				}
			}
		}
	else {
		double complex tempvector[2*(long unsigned) j +1];
		double complex *temp = tempvector+(long unsigned) j;

		for (m2 = -(int)j; m2 <= (int)j; m2++)
			for (m1 = -(int)j; m1 <= (int)j; m1++) {
				d = DLM_(j, (double)m1, (double)m2, alpha, beta, gamma); 
	        	temp[m2] += d * initial[m1];
				}

		for (m2 = -(int)j; m2 <= (int)j; m2++)
			final[m2] = temp[m2];
		}
}

