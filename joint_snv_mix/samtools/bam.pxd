'''
Created on 2012-01-17

@author: Andrew Roth
'''
from libc.stdlib cimport free
from libc.string cimport strdup

from joint_snv_mix.samtools.samtools_clib cimport bam_index_load, bam_index_destroy, bam_parse_region, samopen, \
                                                  samclose, samfile_t, bam_index_t

from joint_snv_mix.samtools.pileup cimport PileupIterator                                                 

cdef class BamFile:
    cdef char * _file_name
    cdef samfile_t * _bam_file
    cdef bam_index_t * _index
    
    cdef list _refs
    cdef dict _tids
    
    cdef samfile_t * get_file_pointer(self)
    cdef bam_index_t * get_index(self)
