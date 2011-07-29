from csamtools cimport Samfile, IteratorColumnRegion

from joint_snv_mix.counters.counter cimport Counter, CounterRefIterator, CounterRow
from joint_snv_mix.counters.ref_iterator cimport CRefIterator
from joint_snv_mix.counters.shared cimport counts_struct, base_counts_struct, column_struct, strcmp
    
cdef class BaseCounter(Counter):
    cdef Samfile _bam_file
    cdef int _min_base_qual
    cdef int _min_map_qual
    
cdef class BaseCounterRow(CounterRow):
    cdef counts_struct _counts
    
    cdef int get_counts(self, char * base)
    cdef base_counts_struct get_base_counts(self, char * base)

cdef class BaseCounterRefIterator(CounterRefIterator):
    cdef int _min_base_qual
    cdef int _min_map_qual
    cdef CRefIterator _ref_iter
        
    cdef parse_current_position(self)
    cdef advance_position(self)
    
    cdef counts_struct _parse_column(self, column_struct column)
