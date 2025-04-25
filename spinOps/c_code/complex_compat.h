#ifndef COMPLEX_COMPAT_H
#define COMPLEX_COMPAT_H

/* If we're on a truly C99‐capable compiler (GCC/Clang)…
   but _not_ MSVC or clang-cl in msvc mode: */
// #if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L \
//     && !defined(_MSC_VER)
#if true

/* use the system’s <complex.h> */
#include <complex.h>
#include <math.h>

typedef float _Complex fcomplex;
typedef double _Complex dcomplex;
typedef long double _Complex lcomplex;

/* Constructors/macro aliases (C99 already provides CMPLXF, CMPLX, CMPLXL) */
/* Ensure we have a unified CMPLX for dcomplex */
#ifndef CMPLX
#define CMPLX(r, i) ((dcomplex)((r) + (i) * I))
#endif

/* Arithmetic helpers so fallback code compiles unmodified */
static inline dcomplex cadd(dcomplex a, dcomplex b) { return a + b; }
static inline dcomplex csub(dcomplex a, dcomplex b) { return a - b; }
static inline dcomplex cmul(dcomplex a, dcomplex b) { return a * b; }
static inline dcomplex cdiv(dcomplex a, dcomplex b) { return a / b; }

#else
/* fallback for Windows/MSVC: */
#include <math.h>

/* simple struct-based complex type */
typedef struct
{
  float re, im;
} fcomplex;

typedef struct
{
  double re, im;
} dcomplex;

typedef struct
{
  long double re, im;
} lcomplex;

/* constructors */
static inline fcomplex CMPLXF(float x, float y) { return (fcomplex){x, y}; }
static inline dcomplex CMPLX(double x, double y) { return (dcomplex){x, y}; }
static inline lcomplex CMPLXL(long double x, long double y) { return (lcomplex){x, y}; }

/* conjugate */
static inline fcomplex conjf(fcomplex z) { return (fcomplex){z.re, -z.im}; }
static inline dcomplex conj(dcomplex z) { return (dcomplex){z.re, -z.im}; }
static inline lcomplex conjl(lcomplex z) { return (lcomplex){z.re, -z.im}; }

/* magnitude, phase… */
static inline double cabs(dcomplex z) { return sqrt(z.re * z.re + z.im * z.im); }
static inline double carg(dcomplex z) { return atan2(z.im, z.re); }

/* float complex */
static inline fcomplex caddf(fcomplex a, fcomplex b) {
    return (fcomplex){ a.re + b.re, a.im + b.im };
}
static inline fcomplex csubf(fcomplex a, fcomplex b) {
    return (fcomplex){ a.re - b.re, a.im - b.im };
}
static inline fcomplex cmulf(fcomplex a, fcomplex b) {
    return (fcomplex){
        a.re * b.re - a.im * b.im,
        a.re * b.im + a.im * b.re
    };
}
static inline fcomplex cdivf(fcomplex a, fcomplex b) {
    float denom = b.re*b.re + b.im*b.im;
    return (fcomplex){
        (a.re*b.re + a.im*b.im)/denom,
        (a.im*b.re - a.re*b.im)/denom
    };
}

static inline dcomplex cadd(dcomplex a, dcomplex b) {
    return (dcomplex){ a.re + b.re, a.im + b.im };
}
static inline dcomplex csub(dcomplex a, dcomplex b) {
    return (dcomplex){ a.re - b.re, a.im - b.im };
}
static inline dcomplex cmul(dcomplex a, dcomplex b) {
    return (dcomplex){ a.re*b.re - a.im*b.im,
                       a.re*b.im + a.im*b.re };
}
static inline dcomplex cdiv(dcomplex a, dcomplex b) {
    double denom = b.re*b.re + b.im*b.im;
    return (dcomplex){ (a.re*b.re + a.im*b.im)/denom,
                       (a.im*b.re - a.re*b.im)/denom };
}

#endif /* C99 vs fallback */

#endif /* COMPLEX_COMPAT_H */
