#include "spin.h"
#include "fact.h"

#define MAX_TWO_I  11    // supports 2I = 1,2,…,11  i.e. I = ½,1,…,11/2
#define MAX_L       8    // supports l = 0,1,…,8

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

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
#include "spin.h"   // brings in small_fac[] and MAX_SMALL_FAC

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
                matrix[bra][ket] = 1. / sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
                matrix[bra][ket] -= 1 / sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
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
                matrix[bra][ket] += 1 / sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
                matrix[bra][ket] -= 1 / sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
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
                matrix[bra][ket] = I / sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
                matrix[bra][ket] += I / sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
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
                matrix[bra][ket] += I / sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
                matrix[bra][ket] += I / sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
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
            matrix[bra][ket] += tlm_(1., 0., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
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
                matrix[bra][ket] += I / sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
                matrix[bra][ket] += I / sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
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
                matrix[bra][ket] = -sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
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
                matrix[bra][ket] = -sqrt(2) * tlm_(1., 1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
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
                matrix[bra][ket] = sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
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
                matrix[bra][ket] = sqrt(2) * tlm_(1., -1., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
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
                matrix[bra][ket] = tlm_(L, M, i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
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
                matrix[bra][ket] = unit_tlm_(L, M, i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}


/*!
 @function get_single_spin_C0_
 */
void get_single_spin_C0_(double complex *operator, int spin_index, int *i_times_2, int total_spin_count)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates = number_of_states_(total_spin_count, i_times_2);
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    float spin = i_times_2[spin_index]/2.0;
    double pi01 = 0.3577708763999664 * (spin*(spin+1) - 0.75);
    double pi03 = 0.848528137423857;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
                matrix[bra][ket] = pi01 * tlm_(1., 0., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
                matrix[bra][ket] += pi03 * tlm_(3., 0., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}


/*!
 @function get_single_spin_C2_
 */
void get_single_spin_C2_(double complex *operator, int spin_index, int *i_times_2, int total_spin_count)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates = number_of_states_(total_spin_count, i_times_2);
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    float spin = i_times_2[spin_index]/2.0;
    double pi21 = 0.1069044967649698 * (spin*(spin+1) - 0.75);
    double pi23 = -1.01418510567422;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
                matrix[bra][ket] = pi21 * tlm_(1., 0., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
                matrix[bra][ket] += pi23 * tlm_(3., 0., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
        }
    }
    free(qnum_data);
}


/*!
 @function get_single_spin_C4_
 */
void get_single_spin_C4_(double complex *operator, int spin_index, int *i_times_2, int total_spin_count)
{
    if (spin_index < 0 || spin_index > total_spin_count - 1)
        return;
    int nstates = number_of_states_(total_spin_count, i_times_2);
    int *qnum_data = createQuantumNumbers(total_spin_count, i_times_2);
    int (*qnum)[nstates] = (int (*)[nstates])qnum_data;
    double complex(*matrix)[nstates] = (double complex(*)[nstates])operator;

    float spin = i_times_2[spin_index]/2.0;
    double pi41 = -0.1434274331201272 * (spin*(spin+1) - 0.75);
    double pi43 = 1.285079208231372;

    for (int bra = 0; bra < nstates; bra++)
    {
        for (int ket = 0; ket < nstates; ket++)
        {
            int del = systemDeltaProduct(qnum_data, total_spin_count, nstates, spin_index, bra, ket);
            if (del == 0)
                matrix[bra][ket] = 0;
            else
                matrix[bra][ket] = pi41 * tlm_(1., 0., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
                matrix[bra][ket] += pi43 * tlm_(3., 0., i_times_2[spin_index], qnum[spin_index][bra], qnum[spin_index][ket]) * del;
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
