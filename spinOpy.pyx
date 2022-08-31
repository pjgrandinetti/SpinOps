from numpy cimport ndarray
import numpy as np

cdef extern from "spinOp.h":
    double Clebsch_(double j1,double m1,double j2,double m2,double j,double m)
    double TLM_(double l,double m,double j1,double m1,double j2,double m2)
    double unitTLM_(double l,double m,double j1,double m1,double j2,double m2)
    int numberOfStates_(int spinCount, int *allSpinsTimesTwo)

def Clebsch(j1: double, m1: double, j2: double, m2: double, j: double, m: double):
    return Clebsch_(j1,m1,j2,m2,j,m)

def TLM(l: double, m: double, j1: double, m1: double, j2: double, m2: double):
    return TLM_(l,m,j1,m1,j2,m2)

def unitTLM(l: double, m: double, j1: double, m1: double, j2: double, m2: double):
    return unitTLM_(l,m,j1,m1,j2,m2)

def numberOfStates(list allSpinsTimesTwo):
    cdef int spinCount = len(allSpinsTimesTwo)
    cdef ndarray[int] spins=np.array(allSpinsTimesTwo,dtype=int)
    return numberOfStates_(spinCount, &spins[0])

