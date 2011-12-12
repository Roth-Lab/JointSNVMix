cdef extern from "stdlib.h":
    int strcmp (char * , char *)
    
cdef extern from "ctype.h":
    int toupper(int c)

#=======================================================================================================================
# Generic counts structs
#=======================================================================================================================
ctypedef struct counts_struct:
    int A
    int C
    int G
    int T

ctypedef struct base_counts_struct:
    char * base
    int counts

#=======================================================================================================================
# Binary data structs
#=======================================================================================================================
ctypedef struct binary_counts_struct:
    int A
    int B
    
ctypedef struct binary_quality_struct:
    int * A
    int * B

ctypedef struct base_map_qualities_struct:
    binary_counts_struct depth
    binary_quality_struct base_quals
    binary_quality_struct map_quals

#=======================================================================================================================
# Pileup Column struct
#=======================================================================================================================
ctypedef struct column_struct:
    char * ref
    int position
    int depth
    
    char * bases
    int * base_quals
    int * map_quals