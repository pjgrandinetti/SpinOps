cdef extern from "spinOp.h":
    double fac(double x)

def pyFac(x: double):
    return fac(x)
