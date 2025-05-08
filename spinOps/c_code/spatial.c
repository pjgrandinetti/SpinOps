#include "spatial.h"
#include "fact.h"


/* General reduced Wigner-d as before */
double _wigner_d_reduced_general(int two_J,
                                 int two_m1,
                                 int two_m2,
                                 double theta)
{
    double J  = 0.5 * two_J;
    double m1 = 0.5 * two_m1;
    double m2 = 0.5 * two_m2;

    double c2 = cos(0.5 * theta);
    double s2 = sin(0.5 * theta);

    /* summation limits */
    int kmin = 0;
    int tmp  = (two_m1 + two_m2)/2;
    if(tmp < 0) kmin = -tmp;
    int kmax1 = (two_J - two_m1) / 2;
    int kmax2 = (two_J - two_m2) / 2;
    int kmax  = kmax1 < kmax2 ? kmax1 : kmax2;

    double sum = 0.0;
    for(int k = kmin; k <= kmax; ++k) {
        int a = k;
        int b = (two_J - two_m1)/2 - k;
        int c = (two_J - two_m2)/2 - k;
        int d = tmp + k;
        if(b<0||c<0||d<0) continue;

        double sign = (k&1)? -1.0 : 1.0;
        double p1   = pow(c2, m1 + m2 + 2.0*k);
        double p2   = pow(s2, 2.0*J - m1 - m2 - 2.0*k);
        double denom= fac_int(a)*fac_int(b)*fac_int(c)*fac_int(d);
        sum += sign * (p1 * p2)/denom;
    }

    int expo = (two_J - two_m2)/2;
    double phase = (expo&1)? -1.0 : 1.0;

    int jm1    = (two_J + two_m1)/2;
    int jm1bar = (two_J - two_m1)/2;
    int jm2    = (two_J + two_m2)/2;
    int jm2bar = (two_J - two_m2)/2;
    double norm = sqrt(
        fac_int(jm1)*fac_int(jm1bar) *
        fac_int(jm2)*fac_int(jm2bar)
    );

    return phase * norm * sum;
}

/* Special-case for J=2 (two_J==4) using your hard-coded expressions: */
double _wigner_d_reduced_J2(int two_m1, int two_m2, double theta) {
    double cx = cos(theta);
    double sx = sin(theta);
    /* use exactly the same constants you inlined before */
    if (two_m1 ==  4) {
        if (two_m2 ==  4) return  0.25 * (1+cx)*(1+cx);
        if (two_m2 ==  2) return -0.5  * sx*(1+cx);
        if (two_m2 ==  0) return  0.6123724355 * sx*sx;
        if (two_m2 == -2) return -0.5  * sx*(1-cx);
        if (two_m2 == -4) return  0.25 * (1-cx)*(1-cx);
    }
    else if (two_m1 == -4) {
        if (two_m2 ==  4) return  0.25 * (1-cx)*(1-cx);
        if (two_m2 ==  2) return  0.5  * sx*(1-cx);
        if (two_m2 ==  0) return  0.6123724355 * sx*sx;
        if (two_m2 == -2) return  0.5  * sx*(1+cx);
        if (two_m2 == -4) return  0.25 * (1+cx)*(1+cx);
    }
    else if (two_m1 ==  2) {
        if (two_m2 ==  4) return  0.5  * sx*(1+cx);
        if (two_m2 ==  2) return  0.5  * (2*cx*cx + cx - 1);
        if (two_m2 ==  0) return -1.224744871 * sx*cx;
        if (two_m2 == -2) return -0.5  * (2*cx*cx - cx - 1);
        if (two_m2 == -4) return -0.5  * sx*(1-cx);
    }
    else if (two_m1 ==  0) {
        if (two_m2 ==  4 || two_m2 == -4) return  0.6123724355 * sx*sx;
        if (two_m2 ==  2) return  1.224744871 * sx*cx;
        if (two_m2 ==  0) return  1.5 * cx*cx - 0.5;
        if (two_m2 == -2) return -1.224744871 * sx*cx;
    }
    else if (two_m1 == -2) {
        if (two_m2 ==  4) return  0.5  * sx*(1-cx);
        if (two_m2 ==  2) return -0.5  * (2*cx*cx - cx - 1);
        if (two_m2 ==  0) return  1.224744871 * sx*cx;
        if (two_m2 == -2) return  0.5  * (2*cx*cx + cx - 1);
        if (two_m2 == -4) return -0.5  * sx*(1+cx);
    }
    /* (should never fall through for valid two_m1,two_m2) */
    return 0.0;
}

/**
 * Top-level entry point: picks the J=2 shortcut if possible,
 * otherwise falls back to the full summation.
 */
double wigner_d_(const int two_J, const int two_m1, const int two_m2, const double theta)
{
    if (two_J == 4) {
        return _wigner_d_reduced_J2(two_m1, two_m2, theta);
    }
    else {
        return _wigner_d_reduced_general(two_J, two_m1, two_m2, theta);
    }
}

/**
 * Computes the Wigner D-matrix element D^l_{m1,m2}(α, β, γ)
 * using the Euler angle decomposition:
 *     D^l_{m1,m2}(α, β, γ) = e^{-i m1 α} d^l_{m1,m2}(β) e^{-i m2 γ}
 *
 * All quantum numbers are passed as 2× their values.
 *
 * @param two_l   Integer: 2 × l
 * @param two_m1  Integer: 2 × m₁
 * @param two_m2  Integer: 2 × m₂
 * @param alpha   Euler angle α (in radians)
 * @param beta    Euler angle β  (in radians)
 * @param gamma   Euler angle γ  (in radians)
 *
 * @return        Complex D-matrix element D^l_{m1,m2}(α, β, γ)
 */
double complex DLM_(const int two_l, const int two_m1, const int two_m2, const double alpha, const double beta, const double gamma)
{
    // Compute the phase factor exp(-i(m1*α + m2*γ))
    double phase_angle = 0.5 * two_m1 * alpha + 0.5 * two_m2 * gamma;
    double d_element = wigner_d_(two_l, two_m1, two_m2, beta);
    return d_element * cexp(-I * phase_angle);
}

/**
 * Applies a rotation defined by Euler angles (α, β, γ) to a spin state vector
 * of angular momentum j, using Wigner D-matrix elements.
 *
 * All quantum numbers are passed as 2× their values (i.e., integer-doubled).
 *
 * @param two_j     Integer: 2 × j
 * @param initial   Input vector of size (two_j + 1)
 * @param alpha     Euler angle α (in radians)
 * @param beta      Euler angle β  (in radians)
 * @param gamma     Euler angle γ  (in radians)
 * @param final     Output vector (must be preallocated, size = two_j + 1)
 */
void Rot_(const int two_j, const double complex *initial, const double alpha, const double beta, const double gamma, double complex *final)
{
    int length = two_j + 1;

    // Identity rotation shortcut
    if (alpha == 0.0 && beta == 0.0 && gamma == 0.0) {
        for (int i = 0; i < length; i++) {
            final[i] = initial[i];
        }
        return;
    }

    // Handle j == 2 (two_j == 4) using hard-coded Wigner d
    if (two_j == 4) {
        for (int two_m2 = -two_j; two_m2 <= two_j; two_m2 += 2) {
            int index2 = (two_j + two_m2) / 2;
            final[index2] = 0.0 + 0.0 * I;
            for (int two_m1 = -two_j; two_m1 <= two_j; two_m1 += 2) {
                int index1 = (two_j + two_m1) / 2;
                double db = wigner_d_(two_j, two_m1, two_m2, beta);
                double pha = 0.5 * (two_m1 * alpha + two_m2 * gamma);
                double complex d = db * cexp(-I * pha);
                final[index2] += d * initial[index1];
            }
        }
    }
    else {
        // General case
        double complex temp[length];
        for (int i = 0; i < length; i++) temp[i] = 0.0 + 0.0 * I;

        for (int two_m2 = -two_j; two_m2 <= two_j; two_m2 += 2) {
            int index2 = (two_j + two_m2) / 2;
            for (int two_m1 = -two_j; two_m1 <= two_j; two_m1 += 2) {
                int index1 = (two_j + two_m1) / 2;
                double complex D = DLM_(two_j, two_m1, two_m2, alpha, beta, gamma);
                temp[index2] += D * initial[index1];
            }
        }

        for (int i = 0; i < length; i++) {
            final[i] = temp[i];
        }
    }
}

/*!
 @function get_rho1_pas_
 */
void get_rho1_pas_(double complex *tensor, const double zeta) {
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
 @function get_rho2_pas_
 */
void get_rho2_pas_(double complex *tensor, const double zeta, const double eta) {
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

