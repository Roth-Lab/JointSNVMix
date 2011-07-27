cdef extern from "stdlib.h":
    int strcmp (char * , char *)

ctypedef struct counts_struct:
    int A
    int C
    int G
    int T

ctypedef struct base_counts_struct:
    char * base
    int counts

ctypedef struct binary_counts_struct:
    int A
    int B