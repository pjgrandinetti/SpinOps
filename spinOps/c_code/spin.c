#include "spin.h"

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define MAX_SMALL_FAC 32
static const double small_fac[MAX_SMALL_FAC + 1] = {
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

#define MAX_LOGFAC 100

static double logfac_table[MAX_LOGFAC + 1];
static int logfac_initialized = 0;

/**
 * Initialize log(factorial) table up to MAX_LOGFAC.
 */
void init_logfac_table()
{
    logfac_table[0] = 0.0;
    for (int n = 1; n <= MAX_LOGFAC; ++n)
        logfac_table[n] = logfac_table[n - 1] + log((double)n);
    logfac_initialized = 1;
}

/**
 * Retrieve log(n!) using precomputed table.
 */
static inline double logfac(int n)
{
    if (!logfac_initialized)
        init_logfac_table();
    if (n < 0 || n > MAX_LOGFAC)
    {
        fprintf(stderr, "Error: logfac(%d) out of bounds.\n", n);
        return NAN;
    }
    return logfac_table[n];
}

double fac(double x) {
    if (x < 0) return 0;
    double result = 1.0;
    for (int i = 2; i <= (int)x; ++i)
        result *= i;
    return result;
}

static inline double fac_int(int n) {
    if (n < 0) return 0.0;
    if (n <= MAX_SMALL_FAC)
        return small_fac[n];
    // fallback for larger n, though this should not occur
    double result = 1.0;
    for (int i = 2; i <= n; ++i)
        result *= i;
    return result;
}

double clebsch_(const int two_J1, const int two_M1,
                const int two_J2, const int two_M2,
                const int two_J,  const int two_M)
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
    for (int k = kmin; k <= kmax; ++k) {
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

/**
 * Computes ⟨j2 m2 | T^{(l)}_m | j1 m1⟩ using Clebsch–Gordan × RME.
 *
 * @param l        Tensor rank (integer)
 * @param m        Tensor component (integer)
 * @param two_j1   2*j1
 * @param two_m1   2*m1
 * @param two_j2   2*j2
 * @param two_m2   2*m2
 */
double tlm_(const int l, const int m, const int two_j1, const int two_m1, const int two_j2, const int two_m2)
{
    // Check domain validity
    if (two_j1 != two_j2 ||
        abs(two_m1) > two_j1 || abs(two_m2) > two_j2 ||
        abs(m) > l || l < 0 || l > two_j1)
        return 0.0;

    // Clebsch–Gordan coefficient: ⟨j2 m2 l m | j1 m1⟩
    double cgc = clebsch_(two_j2, two_m2, 2 * l, 2 * m, two_j1, two_m1);
    if (cgc == 0.0 || isnan(cgc))
        return 0.0;

    // Reduced matrix element
    int j = two_j1; // = two_j2

    double log_rme;
    if (l <= MAX_SMALL_FAC &&
        j + l + 1 <= MAX_LOGFAC &&
        2 * l <= MAX_LOGFAC &&
        j - l >= 0 && j - l <= MAX_LOGFAC)
    {
        log_rme = 2.0 * log(small_fac[l]) +
                  logfac_table[j + l + 1] -
                  l * log(2.0) -
                  logfac_table[2 * l] -
                  logfac_table[j - l];
    }
    else
    {
        log_rme = 2.0 * logfac(l) +
                  logfac(j + l + 1) -
                  l * log(2.0) -
                  logfac(2 * l) -
                  logfac(j - l);
    }

    double rme = exp(0.5 * log_rme);
    return cgc * rme;
}

/**
 * Compute ⟨j2 m2 | T̂_{l,m} | j1 m1⟩
 * using orthonormal unit irreducible spherical tensor operators.
 *
 * All angular momentum arguments (j1, j2, m1, m2) must be passed as 2× their values (integers).
 * Arguments l and m are passed as standard integers.
 */
double unit_tlm_(const int l, const int m, const int two_j1, const int two_m1, const int two_j2, const int two_m2)
{
    if (two_j1 != two_j2 ||
        abs(two_m1) > two_j1 || abs(two_m2) > two_j2 ||
        abs(m) > l || l < 0)
        return 0.0;

    // Clebsch–Gordan coefficient ⟨j2 m2 l m | j1 m1⟩
    double cg = clebsch_(two_j2, two_m2, 2 * l, 2 * m, two_j1, two_m1);

    if (cg == 0.0 || isnan(cg))
        return 0.0;

    // Normalization factor from unit tensor definition
    return cg / sqrt((double)(two_j1 + 1));
}

/**
 * Compute the total number of basis states for a multi-spin system.
 * Each spin has 2I+1 states (I given as 2*I in i_times_2).
 */
int number_of_states_(int total_spin_count, const int *i_times_2)
{
    int nstates = 1;
    for (int i = 0; i < total_spin_count; i++)
    {
        nstates *= i_times_2[i] + 1;
    }
    return nstates;
}

/**
 * Create an array of quantum numbers [2m_i] for each spin index and basis state.
 *
 * @param total_spin_count  Number of spins
 * @param i_times_2         Array of 2*I values for each spin
 * @return Pointer to flattened array: layout is [total_spin_count][nstates]
 *         Each column is a spin configuration: [2m_0, 2m_1, ..., 2m_{N-1}]
 */
int *createQuantumNumbers(int total_spin_count, const int *i_times_2)
{
    const int nstates = number_of_states_(total_spin_count, i_times_2);

    // Allocate matrix [total_spin_count][nstates] in row-major layout
    int *qnum_data = malloc(sizeof(int) * total_spin_count * nstates);
    if (!qnum_data)
    {
        fprintf(stderr, "Error: memory allocation failed in createQuantumNumbers.\n");
        return NULL;
    }

    int *current_state = calloc(total_spin_count, sizeof(int));
    if (!current_state)
    {
        free(qnum_data);
        fprintf(stderr, "Error: memory allocation failed for current_state.\n");
        return NULL;
    }

    for (int s = 0; s < nstates; s++)
    {
        for (int i = 0; i < total_spin_count; i++)
        {
            int two_I = i_times_2[i];
            int n_levels = two_I + 1;
            int two_m = -two_I + 2 * current_state[i];
            qnum_data[i * nstates + s] = two_m;
        }

        // Increment current_state[] like a mixed-base counter
        for (int i = total_spin_count - 1; i >= 0; i--)
        {
            if (++current_state[i] <= i_times_2[i])
                break;
            current_state[i] = 0;
        }
    }

    free(current_state);
    return qnum_data;
}

/**
 * Integer delta function δ_{m1,m2}, where inputs are 2×m values.
 *
 * @param two_m1  Integer: 2 × m₁
 * @param two_m2  Integer: 2 × m₂
 * @return        1 if m₁ == m₂, 0 otherwise
 */
static inline int deltaFunction(const int two_m1, const int two_m2)
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
int systemDeltaProduct(const int *qnum_data,
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

/*!
 @function get_single_spin_Ix_
 */
void get_single_spin_Ix_(double complex *operator, int spin_index, int *i_times_2, int total_spin_count)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates = number_of_states_(total_spin_count, i_times_2);
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
            {
                matrix[bra][ket] = 1 / sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
                matrix[bra][ket] -= 1 / sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

/*!
 @function get_Ix_
 */
void get_Ix_(double complex *operator, int *spin_indexes, int spin_count, int *i_times_2, int total_spin_count)
{
    for (int i = 0; i < spin_count; i++)
        if (spin_indexes[i] < 0 || spin_indexes[i] > total_spin_count - 1)
            return;

    int nstates = number_of_states_(total_spin_count, i_times_2);
    memset(operator, 0, nstates * nstates * sizeof(double complex)); // zero operator
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            for (int i = 0; i < spin_count; i++)
            {
                int spin_index = spin_indexes[i];
                int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
                matrix[bra][ket] += 1 / sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
                matrix[bra][ket] -= 1 / sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

/*!
 @function get_single_spin_Iy_
 */
void get_single_spin_Iy_(double complex *operator, int spin_index, int *i_times_2, int total_spin_count)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates = number_of_states_(total_spin_count, i_times_2);
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
            {
                matrix[bra][ket] = I / sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
                matrix[bra][ket] += I / sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

/*!
 @function get_Iy_
 */
void get_Iy_(double complex *operator, int *spin_indexes, int spin_count, int *i_times_2, int total_spin_count)
{
    for (int i = 0; i < spin_count; i++)
        if (spin_indexes[i] < 0 || spin_indexes[i] > total_spin_count - 1)
            return;

    int nstates = number_of_states_(total_spin_count, i_times_2);
    memset(operator, 0, nstates * nstates * sizeof(double complex)); // zero operator
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            for (int i = 0; i < spin_count; i++)
            {
                int spin_index = spin_indexes[i];
                int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
                matrix[bra][ket] += I / sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
                matrix[bra][ket] += I / sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

/*!
 @function get_single_spin_Iz_
 */
void get_single_spin_Iz_(double complex *operator, int spin_index, int *i_times_2, int total_spin_count)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates = number_of_states_(total_spin_count, i_times_2);
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            matrix[bra][ket] += tlm_(1., 0., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function get_Iz_
 */
void get_Iz_(double complex *operator, int *spin_indexes, int spin_count, int *i_times_2, int total_spin_count)
{
    for (int i = 0; i < spin_count; i++)
        if (spin_indexes[i] < 0 || spin_indexes[i] > total_spin_count - 1)
            return;

    int nstates = number_of_states_(total_spin_count, i_times_2);
    memset(operator, 0, nstates * nstates * sizeof(double complex)); // zero operator
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            for (int i = 0; i < spin_count; i++)
            {
                int spin_index = spin_indexes[i];
                int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
                matrix[bra][ket] += I / sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
                matrix[bra][ket] += I / sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

/*!
 @function get_single_spin_Ip_
 */
void get_single_spin_Ip_(double complex *operator, int spin_index, int *i_times_2, int total_spin_count)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates = number_of_states_(total_spin_count, i_times_2);
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
                matrix[bra][ket] = -sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function get_Ip_
 */
void get_Ip_(double complex *operator, int *spin_indexes, int spin_count, int *i_times_2, int total_spin_count)
{
    for (int i = 0; i < spin_count; i++)
        if (spin_indexes[i] < 0 || spin_indexes[i] > total_spin_count - 1)
            return;

    int nstates = number_of_states_(total_spin_count, i_times_2);
    memset(operator, 0, nstates * nstates * sizeof(double complex)); // zero operator
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            for (int i = 0; i < spin_count; i++)
            {
                int spin_index = spin_indexes[i];
                int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
                matrix[bra][ket] = -sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

/*!
 @function get_single_spin_Im_
 */
void get_single_spin_Im_(double complex *operator, int spin_index, int *i_times_2, int total_spin_count)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates = number_of_states_(total_spin_count, i_times_2);
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
                matrix[bra][ket] = sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function get_Im_
 */
void get_Im_(double complex *operator, int *spin_indexes, int spin_count, int *i_times_2, int total_spin_count)
{
    for (int i = 0; i < spin_count; i++)
        if (spin_indexes[i] < 0 || spin_indexes[i] > total_spin_count - 1)
            return;

    int nstates = number_of_states_(total_spin_count, i_times_2);
    memset(operator, 0, nstates * nstates * sizeof(double complex)); // zero operator
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            for (int i = 0; i < spin_count; i++)
            {
                int spin_index = spin_indexes[i];
                int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
                matrix[bra][ket] = sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
            }
        }
    }
    free(qnum_data);
}

/*!
 @function get_single_spin_Tlm_
 */
void get_single_spin_Tlm_(double complex *operator, int spin_index, int *i_times_2, int total_spin_count, int L, int M)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates = number_of_states_(total_spin_count, i_times_2);
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
                matrix[bra][ket] = tlm_(L, M, i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function get_single_spin_Tlm_unit_
 */
void get_single_spin_Tlm_unit_(double complex *operator, int spin_index, int *i_times_2, int total_spin_count, int L, int M)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates = number_of_states_(total_spin_count, i_times_2);
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
                matrix[bra][ket] = unit_tlm_(L, M, i_times_2[spin_index], qnum[spin_index][bra], i_times_2[spin_index], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}

/*!
 @function getEf_
 */
void getEf_(double complex *operator, int r, int s, int *i_times_2, int total_spin_count)
{
    int nstates = number_of_states_(total_spin_count, i_times_2);
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
 @function getIxf_
 */
void getIxf_(double complex *operator, int r, int s, int *i_times_2, int total_spin_count)
{
    int nstates = number_of_states_(total_spin_count, i_times_2);
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
 @function getIyf_
 */
void getIyf_(double complex *operator, int r, int s, int *i_times_2, int total_spin_count)
{
    int nstates = number_of_states_(total_spin_count, i_times_2);
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
 @function getIzf_
 */
void getIzf_(double complex *operator, int r, int s, int *i_times_2, int total_spin_count)
{
    int nstates = number_of_states_(total_spin_count, i_times_2);
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
 @function getIpf_
 */
void getIpf_(double complex *operator, int r, int s, int *i_times_2, int total_spin_count)
{
    int nstates = number_of_states_(total_spin_count, i_times_2);
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
 @function getImf_
 */
void getImf_(double complex *operator, int r, int s, int *i_times_2, int total_spin_count)
{
    int nstates = number_of_states_(total_spin_count, i_times_2);
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
