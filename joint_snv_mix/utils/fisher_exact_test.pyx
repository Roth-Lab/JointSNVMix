# This code was originally from https://github.com/brentp/fishers_exact_test.

#cython: cdivision=True

cdef PValues _factory(double left, double right, double two):
    cdef PValues instance = PValues.__new__(PValues)
    
    instance.left_tail = left
    instance.right_tail = right
    instance.two_tail = two
    
    return instance

cpdef PValues fisher_exact_test(int a_true, int a_false, int b_true, int b_false):
    # convert the a/b groups to study vs population.
    cdef int k = a_true
    cdef int n = a_false + a_true # total in study
    cdef int C = a_true + b_true
    cdef int G = C + a_false + b_false

    cdef int um = imin2(n, C)
    cdef int lm = imax2(0, n + C - G)
    
    cdef double epsilon = 1e-10
    
    cdef PValues pv

    if um == lm:
        pv = _factory(1.0, 1.0, 1.0)
        return pv

    cdef double cutoff = hypergeometric_probability(k, n, C, G)
    cdef double left_tail = 0, right_tail = 0, two_tail = 0
    cdef int i
    cdef double p

    for i in range(lm, um + 1):
        p = hypergeometric_probability(i, n, C, G)

        if i <= k:
            left_tail += p

        if i >= k:
            right_tail += p

        if p < cutoff + epsilon:
            two_tail += p

    left_tail = left_tail if left_tail < 1.0 else 1.0
    
    right_tail = right_tail if right_tail < 1.0 else 1.0
    
    two_tail = two_tail if two_tail < 1.0 else 1.0
    
    pv = _factory(left_tail, right_tail, two_tail)
    
    return pv