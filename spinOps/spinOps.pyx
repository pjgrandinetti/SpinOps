# cython: language_level=3
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

from numpy cimport ndarray
cimport numpy as cnp
import numpy as np

cdef extern from "complex.h":
    pass

cdef extern from "spin.h":
    double clebsch_(double j1,double m1,double j2,double m2,double j,double m)
    double tlm_(double l,double m,double j1,double m1,double j2,double m2)
    double unit_tlm_(double l,double m,double j1,double m1,double j2,double m2)
    int numberOfStates_(int spinCount, int *spinsTimesTwo)
    void getIx_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
    void getIy_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
    void getIz_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
    void getIp_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
    void getIm_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount)
    void getTlm_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount, int L, int M)
    void getTlm_unit_(double complex *operator, int spinIndex, int *spinsTimesTwo, int spinCount, int L, int M)
    void getEf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);
    void getIxf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);
    void getIyf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);
    void getIzf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);
    void getIpf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);
    void getImf_(double complex *operator, int r, int s, int *spinsTimesTwo, int spinCount);


cpdef clebsch(j1: double, m1: double, j2: double, m2: double, j: double, m: double):
    return clebsch_(j1,m1,j2,m2,j,m)

cpdef tlm(l: double, m: double, j1: double, m1: double, j2: double, m2: double):
    return tlm_(l,m,j1,m1,j2,m2)

cpdef unit_tlm(l: double, m: double, j1: double, m1: double, j2: double, m2: double):
    return unit_tlm_(l,m,j1,m1,j2,m2)

cpdef numberOfStates(list spinsTimesTwo):
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    return numberOfStates_(spinCount, &spins[0])

cpdef createIx(spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIx_(&myOp[0,0], spinIndex,  &spins[0], spinCount)
    return myOp

cpdef createIy(spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIy_(&myOp[0,0], spinIndex,  &spins[0], spinCount)
    return myOp

cpdef createIz(spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIz_(&myOp[0,0], spinIndex,  &spins[0], spinCount)
    return myOp

cpdef createIp(spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIp_(&myOp[0,0], spinIndex,  &spins[0], spinCount)
    return myOp

cpdef createIm(spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIm_(&myOp[0,0], spinIndex,  &spins[0], spinCount)
    return myOp

cpdef createTLM(L,M, spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getTlm_(&myOp[0,0], spinIndex,  &spins[0], spinCount, L, M)
    return myOp

cpdef createTLM_unit(L, M, spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getTlm_unit_(&myOp[0,0], spinIndex,  &spins[0], spinCount, L, M)
    return myOp

cpdef createEf(r, s, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getEf_(&myOp[0,0], r, s,  &spins[0], spinCount)
    return myOp

cpdef createIxf(r, s, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIxf_(&myOp[0,0], r, s,  &spins[0], spinCount)
    return myOp

cpdef createIyf(r, s, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIyf_(&myOp[0,0], r, s,  &spins[0], spinCount)
    return myOp

cpdef createIzf(r, s, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIzf_(&myOp[0,0], r, s,  &spins[0], spinCount)
    return myOp

cpdef createIpf(r, s, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIpf_(&myOp[0,0], r, s,  &spins[0], spinCount)
    return myOp

cpdef createImf(r, s, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getImf_(&myOp[0,0], r, s,  &spins[0], spinCount)
    return myOp


cdef extern from "spatial.h":
    void getrho1_pas_(double complex *tensor, double zeta)
    void getrho2_pas_(double complex *tensor, double zeta, double eta)
    double wigner_d_(double l,double m1,double m2,double beta)
    double complex DLM_(double l,double  m1,double m2, double alpha, double beta, double gamma)
    void Rot_(double j, double complex *initial, double alpha, double beta, double gamma, double complex *final)

cpdef createRho1(double zeta):
    cdef cnp.ndarray[double complex, ndim=1] myOp = np.zeros(3, dtype=np.complex128)
    getrho1_pas_(<double complex *> cnp.PyArray_DATA(myOp), zeta)
    return myOp

cpdef createRho2(double zeta, double eta):
    cdef cnp.ndarray[double complex, ndim=1] myOp = np.zeros(5, dtype=np.complex128)
    getrho2_pas_(<double complex *> cnp.PyArray_DATA(myOp), zeta, eta)
    return myOp

cpdef wigner_d(l: double, m1: double, m2: double, beta: double):
    return wigner_d_(l,m1,m2,beta)

cpdef DLM(l: double, m1: double, m2: double, alpha: double, beta: double, gamma: double):
    return DLM_(l, m1, m2, alpha, beta, gamma)

cpdef Rotate(cnp.ndarray[double complex, ndim=1] initial, double alpha, double beta, double gamma):
    cdef cnp.ndarray[double complex, ndim=1] myOp = np.zeros(len(initial), dtype=np.complex128)
    Rot_((len(initial)-1)/2, <double complex *> cnp.PyArray_DATA(initial), alpha, beta, gamma, <double complex *> cnp.PyArray_DATA(myOp))
    return myOp

