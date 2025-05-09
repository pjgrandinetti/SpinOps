// Single compilation unit: include unified API
#include "angular_momentum.h"

// fac_int and small_fac definition from fact.c
static inline double fac_int(int n)
{
    if (n < 0)
        return 0.0;
    if (n <= MAX_SMALL_FAC)
        return small_fac[n];
    double result = 1.0;
    for (int i = 2; i <= n; ++i)
        return result *= i;
    return result;
}

/* The one—and only—definition of the lookup table */
const double small_fac[MAX_SMALL_FAC + 1] = {
    1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0,
    362880.0, 3628800.0, 39916800.0, 479001600.0, 6227020800.0,
    87178291200.0, 1307674368000.0, 20922789888000.0,
    355687428096000.0, 6402373705728000.0, 121645100408832000.0,
    2432902008176640000.0, 51090942171709440000.0,
    1124000727777607680000.0, 25852016738884976640000.0,
    620448401733239439360000.0, 15511210043330985984000000.0,
    403291461126605635584000000.0, 10888869450418352160768000000.0,
    304888344611713860501504000000.0, 8841761993739701954543616000000.0,
    265252859812191058636308480000000.0,
    8222838654177922817725562880000000.0,
    263130836933693530167218012160000000.0
};

// spatial and rotation implementations follow

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
 @function init_rho1_pas_
 */
void init_rho1_pas_(double complex *tensor, const double zeta) {
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
 @function init_rho2_pas_
 */
void init_rho2_pas_(double complex *tensor, const double zeta, const double eta) {
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

// spin operator implementations follow

double clebsch_(const int two_J1, const int two_M1,
                const int two_J2, const int two_M2,
                const int two_J, const int two_M)
{
    // Selection rules (integer version)
    if (two_M1 + two_M2 != two_M ||
        abs(two_M) > two_J ||
        two_J < abs(two_J1 - two_J2) ||
        two_J > two_J1 + two_J2)
        return 0.0;

    // Precompute constants
    int A = (two_J1 + two_J2 - two_J) / 2;
    int B = (two_J1 - two_M1) / 2;
    int C = (two_J2 + two_M2) / 2;
    int E = (two_J - two_J2 + two_M1) / 2;
    int F = (two_J - two_J1 - two_M2) / 2;

    int kmin = MAX(0, MAX((two_J2 - two_J - two_M1) / 2,
                          (two_J1 - two_J + two_M2) / 2));

    int kmax = MIN(A, MIN(B, C));

    // C1 summation
    double C1 = 0.0;
    for (int k = kmin; k <= kmax; ++k)
    {
        int a = A - k;
        int b = B - k;
        int c = C - k;
        int d = k;
        int e = E + k;
        int f = F + k;

        // All indices are guaranteed non-negative due to k bounds
        double denom = fac_int(a) * fac_int(b) * fac_int(c) *
                       fac_int(d) * fac_int(e) * fac_int(f);

        C1 += ((k % 2 == 0) ? 1.0 : -1.0) / denom;
    }

    // Normalization factors
    double C2 = fac_int((two_J1 + two_J2 - two_J) / 2) *
                fac_int((two_J1 - two_J2 + two_J) / 2) *
                fac_int((-two_J1 + two_J2 + two_J) / 2) *
                (double)(two_J + 1) /
                fac_int((two_J1 + two_J2 + two_J) / 2 + 1);

    double C3 = fac_int((two_J + two_M) / 2) *
                fac_int((two_J - two_M) / 2) *
                fac_int((two_J1 + two_M1) / 2) *
                fac_int((two_J1 - two_M1) / 2) *
                fac_int((two_J2 + two_M2) / 2) *
                fac_int((two_J2 - two_M2) / 2);

    return C1 * sqrt(C2 * C3);
}

/* rme_table[two_I][l] = ⟨I‖T^(l)‖I⟩ for l≤two_I, else 0 */
static double rme_table[MAX_TWO_I+1][MAX_L+1];
/* inv_sqrt2I1[two_I] = 1.0 / √(2I+1) = 1.0 / √(two_I+1) */
static double inv_sqrt2I1[MAX_TWO_I+1];
static int    tables_initialized = 0;

#include <math.h>

static void init_tlm_tables(void)
{
    if (tables_initialized) return;
    for (int two_I = 1; two_I <= MAX_TWO_I; ++two_I) {
        inv_sqrt2I1[two_I] = 1.0 / sqrt((double)(two_I + 1));
        for (int l = 0; l <= MAX_L; ++l) {
            if (l > two_I) {
                rme_table[two_I][l] = 0.0;
            } else {
                /* numerator = l!·l!·(2I+l+1)! */
                double num = small_fac[l] * small_fac[l]
                           * small_fac[two_I + l + 1];
                /* denominator = 2^l·(2l)!·(2I–l)! */
                double den = ldexp(1.0, l)           /* 2^l */
                           * small_fac[2*l]
                           * small_fac[two_I - l];
                rme_table[two_I][l] = sqrt(num/den);
            }
        }
    }
    tables_initialized = 1;
}

double tlm_(const int l,
            const int m,
            const int two_I,
            const int two_m1,
            const int two_m2)
{
    /* 1) ensure we only ever support our pre-chosen range */
    if (two_I < 1 || two_I > MAX_TWO_I ||
        l     < 0 || l     > MAX_L      ||
        l     > two_I      ||
        abs(m) > l                  ||
        two_m2 + 2*m != two_m1)
    {
        return 0.0;
    }
    /* trivial scalar */
    if (l == 0 && m == 0)
        return 1.0;

    /* 2) init tables on first use */
    init_tlm_tables();

    /* 3) get CG ⟨I,m2;l,m|I,m1⟩ */
    double cg = clebsch_( two_I, two_m2,
                          2*l,   2*m,
                          two_I, two_m1 );
    if (cg == 0.0) return 0.0;

    /* 4) Wigner–Eckart assembly */
    return cg
         * inv_sqrt2I1[two_I]
         * rme_table[two_I][l];
}

/*!
 @brief Compute ⟨I,m1| 𝒯_{l,m} |I,m2⟩, the *unit* irreducible spherical tensor element,
        given the non-unit matrix element tlm_().
 @param l        tensor rank ℓ  (0 ≤ ℓ ≤ 2I)
 @param m        component m     (|m| ≤ ℓ)
 @param two_I    2×I
 @param two_m1   2×m1
 @param two_m2   2×m2
 @return         the unit-tensor matrix element ⟨I,m1|𝒯_{ℓ,m}|I,m2⟩
*/
double unit_tlm_(const int l,
                 const int m,
                 const int two_I,
                 const int two_m1,
                 const int two_m2)
{
    /* 1) get the non-unit element (will be zero if selection rules fail) */
    double raw = tlm_(l, m, two_I, two_m1, two_m2);
    if (raw == 0.0)
        return 0.0;

    /* 2) compute the scaling factor:
         scale = (1/ℓ!) * sqrt[ (2ℓ+1)*(2I-ℓ)!*2^ℓ*(2ℓ)! / (2I+ℓ+1)! ]
    */
    double inv_l_fact = 1.0 / fac_int(l);

    /* all factorial arguments here are <= 2I+ℓ+1,
       which for I≤11/2 and ℓ≤8 stays within small_fac[] */
    double num = (2*l + 1)
               * fac_int(two_I - l)
               * ldexp(1.0, l)        /* 2^ℓ */
               * fac_int(2*l);
    double den = fac_int(two_I + l + 1);

    double scale = inv_l_fact * sqrt(num / den);

    /* 3) assemble unit tensor element */
    return raw * scale;
}

/**
 * Compute the total number of basis states for a multi-spin system.
 * Each spin has 2I+1 states (I given as 2*I in two_I).
 */
int number_of_states_(int total_spin_count, const int *two_I)
{
    int nstates = 1;
    for (int i = 0; i < total_spin_count; i++)
    {
        nstates *= two_I[i] + 1;
    }
    return nstates;
}

/**
 * Create an array of quantum numbers for each spin index and basis state.
 *
 * This function allocates and returns a flattened two-dimensional array containing the doubled
 * magnetic quantum numbers (2*m) for each spin and each basis configuration. The array is laid
 * out in row-major order with dimensions [total_spin_count][*nstates], where each column represents
 * one basis state and each row corresponds to a particular spin's 2*m values in that state.
 *
 * Example:
 *   For three spins with I values [1/2, 1, 3/2] (two_I = [1, 2, 3]), if nstates = 12,
 *   the returned array has size 3 × 12, and column 0 contains [−1, −2, −3], column 1 [−1, −2, −1],
 *   and so on, covering all combinations of m values.
 *
 * @param total_spin_count  Number of spins in the system.
 * @param two_I             Input array of length total_spin_count; each entry is 2 × I for that spin.
 * @param[out] nstates      Pointer to an integer that will be set to the total number of basis states
 *                          (the product of (two_I[i] + 1) for i = 0..total_spin_count−1).
 * @return                  Pointer to a malloc'ed array of size total_spin_count × (*nstates).
 *                          Each entry is the 2*m value for a given spin and basis state. Caller
 *                          is responsible for freeing this memory when no longer needed.
 */
int *create_quantum_numbers(int total_spin_count, const int *two_I, int *nstates)
{
    // Compute total number of basis states (product of (2I+1) for each spin)
    *nstates = number_of_states_(total_spin_count, two_I);

    // Allocate flat array [total_spin_count][nstates] in row-major layout
    int *qnum_data = malloc(sizeof(int) * total_spin_count * (*nstates));
    if (!qnum_data)
    {
        fprintf(stderr, "Error: memory allocation failed in create_quantum_numbers.\n");
        return NULL;
    }

    // Mixed-base counter: each element tracks current m index for a spin
    int *current_state = calloc(total_spin_count, sizeof(int));
    if (!current_state)
    {
        free(qnum_data);
        fprintf(stderr, "Error: memory allocation failed for current_state.\n");
        return NULL;
    }

    // Loop over each combined basis state index
    for (int s = 0; s < (*nstates); s++)
    {
        // For each spin, compute doubled m-value: two_m = -two_I[i] + 2 * state_index
        for (int i = 0; i < total_spin_count; i++)
        {
            int two_m = -two_I[i] + 2 * current_state[i];
            // Store at row 'i', column 's' in the flat array
            qnum_data[i * (*nstates) + s] = two_m;
        }

        // Increment mixed-base counter: advance least significant spin index
        for (int i = total_spin_count - 1; i >= 0; i--)
        {
            // Move to next m-value for spin i
            if (++current_state[i] <= two_I[i])
                break;              // no carry needed
            current_state[i] = 0;   // reset and carry to next spin
        }
    }

    // Clean up temporary counter and return the quantum number table
    free(current_state);
    return qnum_data;
}

/**
 * Allocate and initialize a quantum_numbers_t struct.
 * Copies two_I array and builds qnum_data internally.
 */
quantum_numbers_t *create_quantum_numbers_struct(int total_spin_count, const int *two_I)
{
    // Allocate the struct
    quantum_numbers_t *qn = malloc(sizeof(*qn));
    if (!qn) {
        fprintf(stderr, "Error: allocation failed for quantum_numbers_t.\n");
        return NULL;
    }
    // Store spin count
    qn->total_spin_count = total_spin_count;
    // Copy two_I values
    qn->two_I = malloc(total_spin_count * sizeof(int));
    if (!qn->two_I) {
        fprintf(stderr, "Error: allocation failed for two_I copy.\n");
        free(qn);
        return NULL;
    }
    memcpy(qn->two_I, two_I, total_spin_count * sizeof(int));
    // Build quantum number table
    qn->qnum_data = create_quantum_numbers(total_spin_count, qn->two_I, &qn->nstates);
    if (!qn->qnum_data) {
        fprintf(stderr, "Error: create_quantum_numbers failed inside struct init.\n");
        free(qn->two_I);
        free(qn);
        return NULL;
    }
    return qn;
}

/**
 * Free a quantum_numbers_t struct and internal data.
 */
void free_quantum_numbers_struct(quantum_numbers_t *qn)
{
    if (!qn) return;
    free(qn->two_I);
    free(qn->qnum_data);
    free(qn);
}

/**
 * Integer delta function δ_{m1,m2}, where inputs are 2×m values.
 *
 * @param two_m1  Integer: 2 × m₁
 * @param two_m2  Integer: 2 × m₂
 * @return        1 if m₁ == m₂, 0 otherwise
 */
static inline int delta_function(const int two_m1, const int two_m2)
{
    return (two_m1 == two_m2);
}

/**
 * Computes δ_{bra,ket} over all spins except iskip, using 2×m quantum numbers.
 *
 * @param qnum_data         Pointer to [total_spin_count][nstates] matrix of 2×m values
 * @param total_spin_count  Number of spins
 * @param nstates           Number of basis states
 * @param iskip             Spin index to exclude from comparison
 * @param bra               Index of bra state (column index)
 * @param ket               Index of ket state (column index)
 * @return                  1 if states match on all spins except iskip, 0 otherwise
 */
int system_delta_product(const int *qnum_data,
                       const int total_spin_count,
                       const int nstates,
                       const int iskip,
                       const int bra,
                       const int ket)
{
    const int (*qnum)[nstates] = (const int (*)[nstates])qnum_data;
    for (int iSpin = 0; iSpin < total_spin_count; ++iSpin)
    {
        if (iSpin == iskip)
            continue;
        if (qnum[iSpin][bra] != qnum[iSpin][ket])
            return 0;
    }
    return 1;
}

/**
 * Initializes the single‐spin Ix operator for spin_index,
 * using a precomputed quantum_numbers_t.
 *
 * @param op   Preallocated nstates×nstates complex matrix.
 * @param qn   Pointer to quantum_numbers_t (contains two_I, qnum_data, nstates).
 * @param spin_index  Which spin (0 ≤ spin_index < qn->total_spin_count).
 */
void init_single_spin_Ix(double complex *op,
                                int spin_index,
                                const quantum_numbers_t *qn)
{
    int S = qn->total_spin_count;
    int n = qn->nstates;
    const int *two_I = qn->two_I;
    // bounds‐check
    if (spin_index < 0 || spin_index >= S) return;

    // alias for clarity
    const int (*qnum)[n] = (const int (*)[n]) qn->qnum_data;  // VLA for qnum_data
    double complex (*mat)[n] = (double complex (*)[n])op;

    // zero the entire operator
    memset(op, 0, n*n * sizeof *op);

    // loop over basis states
    for (int bra = 0; bra < n; bra++) {
        for (int ket = 0; ket < n; ket++) {
            // only nonzero when all other spins match
            bool other_match = true;
            for (int j = 0; j < S; j++) {
                if (j != spin_index && qnum[j][bra] != qnum[j][ket]) {
                    other_match = false;
                    break;
                }
            }
            if (!other_match)
                continue;

            int m1 = qnum[spin_index][bra];
            int m2 = qnum[spin_index][ket];
            // build Ix = (T_{1,−1} − T_{1,+1})/√2
            double c1 = tlm_(1, -1, two_I[spin_index], m1, m2);
            double c2 = tlm_(1, +1, two_I[spin_index], m1, m2);
            mat[bra][ket] = (c1 - c2) / M_SQRT2;
        }
    }
}

/*! @function init_Ix_ */
void init_Ix_(double complex *operator, int *spin_indexes, int spin_count, const quantum_numbers_t *qn)
{
    int S = qn->total_spin_count;
    int n = qn->nstates;
    const int *two_I = qn->two_I;
    const int (*qnum)[n] = (const int (*)[n])qn->qnum_data;
    double complex (*mat)[n] = (double complex(*)[n])operator;
    // zero the operator
    memset(operator, 0, n*n * sizeof *operator);
    // temporary buffer for single-spin operator
    double complex *tmp = malloc(n*n * sizeof *tmp);
    if (!tmp) return;
    for (int i = 0; i < spin_count; ++i) {
        int spin_index = spin_indexes[i];
        if (spin_index < 0 || spin_index >= S) continue;
        memset(tmp, 0, n*n * sizeof *tmp);
        init_single_spin_Ix(tmp, spin_index, qn);
        for (int bra = 0; bra < n; ++bra) {
            for (int ket = 0; ket < n; ++ket) {
                mat[bra][ket] += ((double complex (*)[n])tmp)[bra][ket];
            }
        }
    }
    free(tmp);
}

/*! @function init_single_spin_Iy */
void init_single_spin_Iy(double complex *operator, int spin_index, const quantum_numbers_t *qn)
{
    int S = qn->total_spin_count;
    int n = qn->nstates;
    const int *two_I = qn->two_I;
    
    // bounds-check
    if (spin_index < 0 || spin_index >= S) return;

    const int (*qnum)[n] = (const int (*)[n])qn->qnum_data;
    double complex (*mat)[n] = (double complex (*)[n])operator;

    // zero the entire operator
    memset(operator, 0, n * n * sizeof *operator);

    // loop over basis states
    for (int bra = 0; bra < n; bra++) {
        for (int ket = 0; ket < n; ket++) {
            // only nonzero when all other spins match (like in init_single_spin_Ix)
            bool other_match = true;
            for (int j = 0; j < S; j++) {
                if (j != spin_index && qnum[j][bra] != qnum[j][ket]) {
                    other_match = false;
                    break;
                }
            }
            if (!other_match)
                continue;

            int m1 = qnum[spin_index][bra];
            int m2 = qnum[spin_index][ket];
            
            // build Iy = i * (T_{1,1} + T_{1,-1}) / sqrt(2)
            // T_{1,-1} corresponds to tlm_(1, -1, ...)
            // T_{1,1}  corresponds to tlm_(1,  1, ...)
            double tlm_minus1 = tlm_(1, -1, two_I[spin_index], m1, m2);
            double tlm_plus1  = tlm_(1,  1, two_I[spin_index], m1, m2);
            
            mat[bra][ket] = I * (tlm_plus1 + tlm_minus1) / M_SQRT2;
        }
    }
}

/*! @function init_Iy_ */
void init_Iy_(double complex *operator, int *spin_indexes, int spin_count, const quantum_numbers_t *qn)
{
    int S = qn->total_spin_count;
    int n = qn->nstates;
    double complex (*mat)[n] = (double complex(*)[n])operator;
    // zero operator
    memset(operator, 0, n*n * sizeof *operator);
    // buffer for single-spin
    double complex *tmp = malloc(n*n * sizeof *tmp);
    if (!tmp) return;
    for (int i = 0; i < spin_count; ++i) {
        int spin_index = spin_indexes[i];
        if (spin_index < 0 || spin_index >= S) continue;
        memset(tmp, 0, n*n * sizeof *tmp);
        init_single_spin_Iy(tmp, spin_index, qn);
        for (int bra = 0; bra < n; ++bra) {
            for (int ket = 0; ket < n; ++ket) {
                mat[bra][ket] += ((double complex(*)[n])tmp)[bra][ket];
            }
        }
    }
    free(tmp);
}

/*!
 @function init_single_spin_Iz_
 */
void init_single_spin_Iz_(double complex *operator, int spin_index, int *two_I, int total_spin_count)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates;
    int *qnum_data = create_quantum_numbers(total_spin_count, two_I, &nstates);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = system_delta_product(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            matrix[bra][ket] += tlm_(1., 0., two_I[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function init_Iz_
 */
void init_Iz_(double complex *operator, int *spin_indexes, int spin_count, int *two_I, int total_spin_count)
{
    for (int i = 0; i < spin_count; i++)
        if (spin_indexes[i] < 0 || spin_indexes[i] > total_spin_count - 1)
            return;

    int nstates;
    int *qnum_data = create_quantum_numbers(total_spin_count, two_I, &nstates);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;

    memset(operator, 0, nstates * nstates * sizeof(double complex)); // zero operator
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            for (int i = 0; i < spin_count; i++)
            {
                int spin_index = spin_indexes[i];
                int del = system_delta_product(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
                matrix[bra][ket] += I / sqrt(2) * tlm_(1., -1., two_I[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
                matrix[bra][ket] += I / sqrt(2) * tlm_(1., 1., two_I[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

/*!
 @function init_single_spin_Ip_
 */
void init_single_spin_Ip_(double complex *operator, int spin_index, int *two_I, int total_spin_count)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates;
    int *qnum_data = create_quantum_numbers(total_spin_count, two_I, &nstates);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = system_delta_product(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
                matrix[bra][ket] = -sqrt(2) * tlm_(1., 1., two_I[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function init_Ip_
 */
void init_Ip_(double complex *operator, int *spin_indexes, int spin_count, int *two_I, int total_spin_count)
{
    for (int i = 0; i < spin_count; i++)
        if (spin_indexes[i] < 0 || spin_indexes[i] > total_spin_count - 1)
            return;

    int nstates;
    int *qnum_data = create_quantum_numbers(total_spin_count, two_I, &nstates);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;

    memset(operator, 0, nstates * nstates * sizeof(double complex)); // zero operator
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            for (int i = 0; i < spin_count; i++)
            {
                int spin_index = spin_indexes[i];
                int del = system_delta_product(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
                matrix[bra][ket] = -sqrt(2) * tlm_(1., 1., two_I[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

/*!
 @function init_single_spin_Im_
 */
void init_single_spin_Im_(double complex *operator, int spin_index, int *two_I, int total_spin_count)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates;
    int *qnum_data = create_quantum_numbers(total_spin_count, two_I, &nstates);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = system_delta_product(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
                matrix[bra][ket] = sqrt(2) * tlm_(1., -1., two_I[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function init_Im_
 */
void init_Im_(double complex *operator, int *spin_indexes, int spin_count, int *two_I, int total_spin_count)
{
    for (int i = 0; i < spin_count; i++)
        if (spin_indexes[i] < 0 || spin_indexes[i] > total_spin_count - 1)
            return;

    int nstates;
    int *qnum_data = create_quantum_numbers(total_spin_count, two_I, &nstates);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;

    memset(operator, 0, nstates * nstates * sizeof(double complex)); // zero operator
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            for (int i = 0; i < spin_count; i++)
            {
                int spin_index = spin_indexes[i];
                int del = system_delta_product(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
                matrix[bra][ket] = sqrt(2) * tlm_(1., -1., two_I[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

/*!
 @function init_single_spin_Tlm_
 */
void init_single_spin_Tlm_(double complex *operator, int spin_index, int *two_I, int total_spin_count, int L, int M)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates;
    int *qnum_data = create_quantum_numbers(total_spin_count, two_I, &nstates);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = system_delta_product(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
                matrix[bra][ket] = tlm_(L, M, two_I[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function init_single_spin_Tlm_unit_
 */
void init_single_spin_Tlm_unit_(double complex *operator, int spin_index, int *two_I, int total_spin_count, int L, int M)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates;
    int *qnum_data = create_quantum_numbers(total_spin_count, two_I, &nstates);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = system_delta_product(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
                matrix[bra][ket] = unit_tlm_(L, M, two_I[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/* generic C_n computation for single spin */
static void init_single_spin_C_generic(double complex *operator, int spin_index, int *two_I, int total_spin_count, double base1, double base3)
{
    if (spin_index < 0 || spin_index >= total_spin_count)
        return;
    int nstates;
    int *qnum_data = create_quantum_numbers(total_spin_count, two_I, &nstates);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex (*matrix)[nstates] = (double complex (*)[nstates])operator;

    double spin = two_I[spin_index] / 2.0;
    double pi1 = base1 * (spin * (spin + 1) - 0.75);
    double pi3 = base3;

    for (int bra = 0; bra < nstates; ++bra) {
        for (int ket = 0; ket < nstates; ++ket) {
            int del = system_delta_product(qnum_data, total_spin_count, nstates,
                                         spin_index, bra, ket);
            if (!del) {
                matrix[bra][ket] = 0;
            } else {
                matrix[bra][ket] = pi1 * tlm_(1., 0., two_I[spin_index],
                                            qnum[spin_index][bra], qnum[spin_index][ket]) * del;
                matrix[bra][ket] += pi3 * tlm_(3., 0., two_I[spin_index],
                                             qnum[spin_index][bra], qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

void init_single_spin_C0_(double complex *operator, int spin_index, int *two_I, int total_spin_count)
{
    init_single_spin_C_generic(operator, spin_index, two_I, total_spin_count,
                                0.3577708763999664, 0.848528137423857);
}

void init_single_spin_C2_(double complex *operator, int spin_index, int *two_I, int total_spin_count)
{
    init_single_spin_C_generic(operator, spin_index, two_I, total_spin_count,
                                0.1069044967649698, -1.01418510567422);
}

void init_single_spin_C4_(double complex *operator, int spin_index, int *two_I, int total_spin_count)
{
    init_single_spin_C_generic(operator, spin_index, two_I, total_spin_count,
                               -0.1434274331201272, 1.285079208231372);
}

/*!
 @function init_Ef_
 */
void init_Ef_(double complex *operator, int r, int s, int *two_I, int total_spin_count)
{
    int nstates = number_of_states_(total_spin_count, two_I);
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            matrix[bra][ket] = 0;
            if (bra == ket && ket == s)
                matrix[bra][ket] = 1;
            else if (bra == ket && ket == r)
                matrix[bra][ket] = 1;
        }
    }
}

/*!
 @function init_Ixf_
 */
void init_Ixf_(double complex *operator, int r, int s, int *two_I, int total_spin_count)
{
    int nstates = number_of_states_(total_spin_count, two_I);
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            matrix[bra][ket] = 0;
            if ((bra == r) && (ket == s))
                matrix[bra][ket] = .5;
            else if ((bra == s) && (ket == r))
                matrix[bra][ket] = .5;
        }
    }
}

/*!
 @function init_Iyf_
 */
void init_Iyf_(double complex *operator, int r, int s, int *two_I, int total_spin_count)
{
    int nstates = number_of_states_(total_spin_count, two_I);
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            matrix[bra][ket] = 0;
            if (bra == r && ket == s)
                matrix[bra][ket] = .5 * I;
            else if (bra == s && ket == r)
                matrix[bra][ket] = -.5 * I;
        }
    }
}

/*!
 @function init_Izf_
 */
void init_Izf_(double complex *operator, int r, int s, int *two_I, int total_spin_count)
{
    int nstates = number_of_states_(total_spin_count, two_I);
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            matrix[bra][ket] = 0;
            if (bra == ket && ket == s)
                matrix[bra][ket] = .5;
            else if (bra == ket && ket == r)
                matrix[bra][ket] = -.5;
        }
    }
}

/*!
 @function init_Ipf_
 */
void init_Ipf_(double complex *operator, int r, int s, int *two_I, int total_spin_count)
{
    int nstates = number_of_states_(total_spin_count, two_I);
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            matrix[bra][ket] = 0;
            if ((ket == r) && (bra == s))
                matrix[bra][ket] = 1;
        }
    }
}

/*!
 @function init_Imf_
 */
void init_Imf_(double complex *operator, int r, int s, int *two_I, int total_spin_count)
{
    int nstates = number_of_states_(total_spin_count, two_I);
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            matrix[bra][ket] = 0;
            if ((bra == r) && (ket == s))
                matrix[bra][ket] = 1;
        }
    }
}
