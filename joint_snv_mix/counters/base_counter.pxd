from csamtools cimport Samfile, IteratorColumnRegion

from joint_snv_mix.counters.counter cimport Counter
from joint_snv_mix.counters.counter_row cimport SingleSampleCounterRow
from joint_snv_mix.counters.ref_iterator cimport CRefIterator, RefIterator
from joint_snv_mix.counters.shared cimport counts_struct, base_counts_struct, column_struct, strcmp
    
cdef class BaseCounter(Counter):
    cdef Samfile _bam_file
    cdef int _min_base_qual
    cdef int _min_map_qual
    
cdef class BaseCounterRow(SingleSampleCounterRow):
    cdef counts_struct _counts

cdef class BaseCounterRefIterator(RefIterator):
    cdef int _min_base_qual
    cdef int _min_map_qual    
    cdef CRefIterator _ref_iter

    cdef counts_struct _parse_column(self, column_struct column)
