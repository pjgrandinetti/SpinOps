from numpy cimport ndarray
import numpy as np

cdef extern from "spinOps.h":
    double Clebsch_(double j1,double m1,double j2,double m2,double j,double m)
    double TLM_(double l,double m,double j1,double m1,double j2,double m2)
    double unitTLM_(double l,double m,double j1,double m1,double j2,double m2)
    int numberOfStates_(int spinCount, int *spinsTimesTwo)

def Clebsch(j1: double, m1: double, j2: double, m2: double, j: double, m: double):
    return Clebsch_(j1,m1,j2,m2,j,m)

def TLM(l: double, m: double, j1: double, m1: double, j2: double, m2: double):
    return TLM_(l,m,j1,m1,j2,m2)

def unitTLM(l: double, m: double, j1: double, m1: double, j2: double, m2: double):
    return unitTLM_(l,m,j1,m1,j2,m2)

def numberOfStates(list spinsTimesTwo):
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    return numberOfStates_(spinCount, &spins[0])

def createIx(spinIndex, list spinsTimesTwo):
    nstates = numberOfStates(spinsTimesTwo)
    cdef int spinCount = len(spinsTimesTwo)
    cdef ndarray[int] spins=np.array(spinsTimesTwo,dtype=np.int32)
    cdef ndarray[double, ndim=2] operator=np.zeros((nstates,nstates), dtype=numpy.cdouble)
    getIx_(&operator[0], &spinIndex[0],  &spins[0], spinCount)
    return operator
