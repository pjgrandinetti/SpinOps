from numpy cimport ndarray
import numpy as np

cdef extern from "complex.h":
    pass

cdef extern from "spinOps.h":
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

def clebsch(j1: double, m1: double, j2: double, m2: double, j: double, m: double):
    return clebsch_(j1,m1,j2,m2,j,m)

def tlm(l: double, m: double, j1: double, m1: double, j2: double, m2: double):
    return tlm_(l,m,j1,m1,j2,m2)

def unit_tlm(l: double, m: double, j1: double, m1: double, j2: double, m2: double):
    return unit_tlm_(l,m,j1,m1,j2,m2)

def numberOfStates(list spinsTimesTwo):
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    return numberOfStates_(spinCount, &spins[0])

def createIx(spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIx_(&myOp[0,0], spinIndex,  &spins[0], spinCount)
    return myOp

def createIy(spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIy_(&myOp[0,0], spinIndex,  &spins[0], spinCount)
    return myOp

def createIz(spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIz_(&myOp[0,0], spinIndex,  &spins[0], spinCount)
    return myOp

def createIp(spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIp_(&myOp[0,0], spinIndex,  &spins[0], spinCount)
    return myOp

def createIm(spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getIm_(&myOp[0,0], spinIndex,  &spins[0], spinCount)
    return myOp

def createTLM(L,M, spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getTlm_(&myOp[0,0], spinIndex,  &spins[0], spinCount, L, M)
    return myOp

def createTLM_unit(L, M, spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double complex, ndim=2] myOp = np.zeros((nstates,nstates), dtype=np.complex128)
    getTlm_unit_(&myOp[0,0], spinIndex,  &spins[0], spinCount, L, M)
    return myOp

cdef extern from "spatialOps.h":
    void getrho1_pas_(double complex *tensor, double zeta, double eta)
    void getrho2_pas_(double complex *tensor, double zeta, double eta)
    double wigner_d_(double l,double m1,double m2,double beta)
    double complex DLM_(double l,double  m1,double m2, double alpha, double beta, double gamma)

def createRho1(zeta, eta):
    cdef ndarray[double complex, ndim=1] myOp = np.zeros(3, dtype=np.complex128)
    getrho1_pas_(&myOp[0], zeta, eta)
    return myOp

def createRho2(zeta, eta):
    cdef ndarray[double complex, ndim=1] myOp = np.zeros(5, dtype=np.complex128)
    getrho2_pas_(&myOp[0], zeta, eta)
    return myOp

def wigner_d(l: double, m1: double, m2: double, beta: double):
    return wigner_d_(l,m1,m2,beta)

def DLM(l: double, m1: double, m2: double, alpha: double, beta: double, gamma: double):
    return DLM_(l, m1, m2, alpha, beta, gamma)


