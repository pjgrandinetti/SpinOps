cdef extern from "spinOp.h":
    double Clebsch(double j1,double m1,double j2,double m2,double j,double m)
    double TLM(double l,double m,double j1,double m1,double j2,double m2)
    double unitTLM(double l,double m,double j1,double m1,double j2,double m2)

def pyClebsch(j1: double, m1: double, j2: double, m2: double, j: double, m: double):
    return Clebsch(j1,m1,j2,m2,j,m)

def pyTLM(l: double, m: double, j1: double, m1: double, j2: double, m2: double):
    return TLM(l,m,j1,m1,j2,m2)

def pyunitTLM(l: double, m: double, j1: double, m1: double, j2: double, m2: double):
    return unitTLM(l,m,j1,m1,j2,m2)
