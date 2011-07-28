# This code was originally from https://github.com/brentp/fishers_exact_test.

#cython: cdivision=True

# Logarithm of n! with algorithmic approximation
# Reference:
#   Lanczos, C. 'A precision approximation of the gamma function',
#   J. SIAM Numer. Anal., B, 1, 86-96, 1964."
#   http://www.matforsk.no/ola/fisher.htm 
cdef float lngamma(int z):
    cdef double x = 0.1659470187408462e-06 / (z + 7)
    x += 0.9934937113930748e-05 / (z + 6)
    x -= 0.1385710331296526 / (z + 5)
    x += 12.50734324009056 / (z + 4)
    x -= 176.6150291498386 / (z + 3)
    x += 771.3234287757674 / (z + 2)
    x -= 1259.139216722289 / (z + 1)
    x += 676.5203681218835 / (z)
    x += 0.9999999999995183
    return log(x) - 5.58106146679532777 - z + (z - 0.5) * log(z + 6.5)

cdef float lnfactorial(int n):
    return 0 if n < 1 else lngamma(n + 1)

# Logarithm of the number of combinations of 'n' objects taken 'p' at a time
cdef float lncombination (int n, int p):
    return lnfactorial(n) - \
           lnfactorial(p) - \
           lnfactorial(n - p)

cdef float hypergeometric_probability(int i, int n, int C, int G):
    return exp(
      lncombination(C, i) + 
      lncombination(G - C, n - i) - 
      lncombination(G, n)
     )
