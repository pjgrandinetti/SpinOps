cdef extern from "spinOp.h":
    double Clebsch_(double j1,double m1,double j2,double m2,double j,double m)
    double TLM_(double l,double m,double j1,double m1,double j2,double m2)
    double unitTLM_(double l,double m,double j1,double m1,double j2,double m2)

def Clebsch(j1: double, m1: double, j2: double, m2: double, j: double, m: double):
    return Clebsch_(j1,m1,j2,m2,j,m)

def TLM(l: double, m: double, j1: double, m1: double, j2: double, m2: double):
    return TLM_(l,m,j1,m1,j2,m2)

def unitTLM(l: double, m: double, j1: double, m1: double, j2: double, m2: double):
    return unitTLM_(l,m,j1,m1,j2,m2)
