#include "spatial.h"
#include "spin.h"

/*!
 @function getrho1_pas_
 */
void getrho1_pas_(double complex *tensor, const double zeta) {
    if (tensor == NULL) {
        fprintf(stderr, "Error: tensor pointer is NULL.\n");
        return;
    }

    const double SQRT_2 = 1.41421356237;
    tensor[0] = 0;
    tensor[1] = -I * SQRT_2 * zeta;
    tensor[2] = 0;
}

/*!
 @function getrho2_pas_
 */
void getrho2_pas_(double complex *tensor, const double zeta, const double eta) {
    // Validate input
    if (tensor == NULL) {
        fprintf(stderr, "Error: tensor pointer is NULL.\n");
        return;
    }

    // Define constants for clarity
    const double SQRT_6_OVER_2 = 1.224744871391589;  // sqrt(6)/2

    // Initialize the tensor components
    tensor[0] = eta * zeta / 2;               // m = -2
    tensor[1] = 0;                            // m = -1
    tensor[2] = SQRT_6_OVER_2 * zeta;         // m = 0
    tensor[3] = 0;                            // m = +1
    tensor[4] = eta * zeta / 2;               // m = +2
}

/*!
 @function wigner_d_
 */
double wigner_d_(const double l,const double m1,const double m2,const double beta)
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
		double sx = sin(beta/2.);
		double cx = cos(beta/2.);
		double sum = 0.;
		int sign = 1;
	                          
		for (int k = 0; k <= l - m2; k++) {
			double k1 = (int)(l - m2 - k);
			double k2 = (int)(l + m1 - k);
			double k3 = (int)(k + m2 - m1);

			if ( k1 >= 0 && k2 >= 0 && k3 >= 0) {   
			 	int n1 = (int)(2 * l + m1 - m2 - 2 * k);
				int n2 = (int)(m2 - m1 + 2 * k);
				double x = mypow(cx, n1);
				double y = mypow(sx, n2);
				sum += sign * x * y / (fac((double)k1) * fac((double)k2) * fac((double)k3) * fac((double)k)); 
				}
		sign = -sign;
		}
		double f = fac(l+m1) * fac(l-m1) * fac(l+m2) * fac(l-m2);
		f = sqrt(f);
		return(sum * f);
		}
	return(0);
}

/*!
 @function DLM_
 */
double complex DLM_(const double l, const double m1, const double m2, const double alpha, const double beta, const double gamma)
{
    double pha = m1 * alpha + m2 * gamma;
    double db = wigner_d_(l, m1, m2, beta);                                   
    return cos(pha) * db - I * sin(pha) * db;
}

/*!
 @function Rot_
 */
void Rot_(const double j, double complex *initial, const double alpha, const double beta, const double gamma, double complex *final)
{
	double m1, m2;
	double complex d;
	int length = 2*(int)j+1;

	if((alpha==0.)&&(beta==0.)&&(gamma==0.)) {
		for (int index = 0; index < length; index++) {
			final[index] = initial[index];
		}
		return;
		} 
	
	if(j==2) {
		double pha, db;
		for (m2 = -2; m2 <= 2; m2++) {
			int index2 = (int)j + m2;
			final[index2] = 0.;
			for (m1 = -2; m1 <= 2; m1++) {
				int index1 = (int)j + m1;
				db = wigner_d_(j, m1, m2, beta);                                   
				pha = m1 * alpha + m2 * gamma;
				d = cos(pha) * db - I* sin(pha) * db;
	        	final[index2] += d * initial[index1];
				}
			}
		}
	else {
		double complex tempvector[2*(long unsigned) j +1];
		double complex *temp = tempvector+(long unsigned) j;

		for (m2 = -j; m2 <= j; m2 = m2+1.) {
			int index2 = (int)j + m2;
			for (m1 = -j; m1 <= j; m1=m1+1.) {
				int index1 = (int)j + m1;
				d = DLM_(j,m1, m2, alpha, beta, gamma); 
	        	temp[index2] += d * initial[index1];
				}
			}

		for (m2 = -j; m2 <= j; m2=m2+1) {
			int index2 = (int)j + m2;
			final[index2] = temp[index2];
		}
	}
}

