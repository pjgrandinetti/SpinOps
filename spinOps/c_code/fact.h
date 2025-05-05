/* fact.h */
#ifndef FACT_H
#define FACT_H

/* Maximum n for which we provide a precomputed table lookup */
#define MAX_SMALL_FAC 32

/* 
 * The lookup table is defined once in fact.c. 
 * Every TU sees it as an extern.
 */
extern const double small_fac[MAX_SMALL_FAC + 1];

/**
 * Compute n! for integer n.
 *  - If n < 0: returns 0.0.
 *  - If 0 ≤ n ≤ MAX_SMALL_FAC: returns small_fac[n].
 *  - If n > MAX_SMALL_FAC: computes by simple multiplication.
 *
 * Defined static inline so that it may be inlined into each TU.
 */
static inline double fac_int(int n)
{
    if (n < 0)
        return 0.0;
    if (n <= MAX_SMALL_FAC)
        return small_fac[n];
    double result = 1.0;
    for (int i = 2; i <= n; ++i)
        result *= i;
    return result;
}

#endif /* FACT_H */
