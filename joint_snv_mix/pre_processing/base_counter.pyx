'''
Created on 2011-06-21

@author: Andrew Roth
'''
from csamtools cimport Samfile, IteratorColumnRegion, PileupProxy, bam1_t, bam1_seq, bam1_qual, bam_pileup1_t, bam_dup1, bam_destroy1

cdef extern from "stdlib.h":
    int strcmp (char * , char *)
    
cdef extern from "stdint.h":
    ctypedef int uint8_t

DEF ASCII_OFFSET = 33

ctypedef struct counts_struct:
    int A
    int C
    int G
    int T

cdef class BaseCounter(object):
    cdef Samfile _bam_file
    cdef int min_base_qual
    cdef int min_map_qual

    def __init__(self, Samfile bam_file, int min_base_qual=10, int min_map_qual=10):
        self._bam_file = bam_file
        
        self.min_base_qual = min_base_qual
        self.min_map_qual = min_map_qual
    
    def iter_ref(self, ref):
        iter = BaseCounterIter(ref, self._bam_file.pileup(ref), self.min_base_qual, self.min_map_qual)
        
        return iter
    
    property refs:
        def __get__(self):
            return self._bam_file.references

cdef class BaseCounterIter(object):
    cdef char * ref
    cdef IteratorColumnRegion pileup_iter
    cdef int min_base_qual
    cdef int min_map_qual

    def __init__(self, char * ref, IteratorColumnRegion pileup_iter, int min_base_qual, int min_map_qual):
        self.ref = ref
        self.pileup_iter = pileup_iter
        self.min_base_qual = min_base_qual
        self.min_map_qual = min_map_qual
    
    def __iter__(self):
        return self
    
    def __next__(self):
        cdef counts_struct counts
        cdef int position
        cdef PileupProxy pileup_column
        cdef BaseCounterRow row
        
        pileup_column = self.pileup_iter.next()
    
        counts = self.parse_pileup_column(pileup_column)
    
        position = pileup_column.pos + 1
    
        row = makeBaseCounterRow(self.ref, position, counts)
        
        return row
    
    cdef counts_struct parse_pileup_column(self, PileupProxy pileup_column):
        cdef int x, qpos, map_qual, base_qual
        cdef char * base
        cdef bam1_t * alignment
        cdef bam_pileup1_t * pileup
        cdef counts_struct counts
        
        counts.A = 0
        counts.C = 0
        counts.G = 0
        counts.T = 0
        
        for x from 0 <= x < pileup_column.n_pu:
            pileup = & pileup_column.plp[x]
            
            if pileup.is_del:
                continue
            
            qpos = pileup.qpos
            
            alignment = bam_dup1(pileup.b) 
            map_qual = alignment.core.qual
            
            if map_qual < self.min_map_qual:
                bam_destroy1(alignment)
                continue            
            
            base_qual = get_qual(alignment, qpos)            
            
            if base_qual < self.min_base_qual:
                bam_destroy1(alignment)
                continue
            
            base = get_base(alignment, qpos)
            
            if strcmp(base, "A") == 0:
                counts.A += 1        
            elif strcmp(base, "C") == 0:
                counts.C += 1
            elif strcmp(base, "G") == 0:
                counts.G += 1
            elif strcmp(base, "T") == 0:
                counts.T += 1
            
            bam_destroy1(alignment)
            
        return counts
    
cdef class BaseCounterRow
cdef BaseCounterRow makeBaseCounterRow(char * ref, int position, counts_struct counts):
     cdef BaseCounterRow row = BaseCounterRow.__new__(BaseCounterRow)
    
     row._ref = ref
     row._position = position
     row._counts = counts
     
     return row

cdef class BaseCounterRow(object):
    cdef char * _ref
    cdef int _position
    cdef counts_struct _counts
    
    def __init__(self):
        raise TypeError("This class cannot be instantiated from Python")
    
    def __str__(self):
        return "\t".join((
                         self._ref,
                         str(self._position),
                         str(self._counts.A),
                         str(self._counts.C),
                         str(self._counts.G),
                         str(self._counts.T)
                         ))
    
    property counts:
        def __get__(self):
            return {
                    'A' : self._counts.A,
                    'C' : self._counts.C,
                    'G' : self._counts.G,
                    'T' : self._counts.T
                    }
    
    property ref:
        def __get__(self):
            return self._ref
    
    property position:
        def __get__(self):
            return self._position
    
    property depth:
        def __get__(self):
            cdef int depth
            
            depth = self._counts.A + self._counts.C + self._counts.G + self._counts.T
            
            return depth

#===============================================================================
# Modified pysam code
#===============================================================================
cdef char * bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"

cdef char * get_base(bam1_t * src, int pos):
    cdef uint8_t * p
    cdef char base
    cdef char[2] base_str

    if src.core.l_qseq == 0: 
        return None
    
    if not src.core.l_qseq:
        return None

    seq = bam1_seq(src)
    
    base = bam_nt16_rev_table[seq[pos / 2] >> 4 * (1 - pos % 2) & 0xf]
    
    base_str[0] = base
    base_str[1] = '\0'
    
    return base_str

cdef int get_qual(bam1_t * src, int pos):    
    cdef uint8_t * p

    p = bam1_qual(src)
    if p[0] == 0xff:
        return None

    return p[pos]
