cdef extern from "math.h":
    double log(double)
    double exp(double)

cdef float lngamma(int z)
cdef float lnfactorial(int n)
cdef float lncombination(int n, int p)
cdef float hypergeometric_probability(int i, int n, int C, int G)
