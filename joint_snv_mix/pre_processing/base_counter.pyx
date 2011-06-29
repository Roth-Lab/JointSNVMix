'''
Use the convention every position at C-level is 0-based. If the position is
accessed by property mechanism from python layer it should be one based.

Created on 2011-06-21

@author: Andrew Roth
'''
from csamtools cimport Samfile, Fastafile, IteratorColumnRegion, PileupProxy, bam1_t, bam1_seq, bam1_qual, bam_pileup1_t, bam_dup1, bam_destroy1

cdef extern from * :
    ctypedef char const_char "const char"
    ctypedef void const_void "const void"

cdef extern from "stdlib.h":
    int strcmp (char * , char *)
    void qsort (void * ARRAY, size_t COUNT, size_t SIZE, int (*COMPARE)(const_void * , const_void *))
    
cdef extern from "stdint.h":
    ctypedef int uint8_t

DEF ASCII_OFFSET = 33

ctypedef struct counts_struct:
    int A
    int C
    int G
    int T

ctypedef struct binary_counts_struct:
    int A
    int B
    
ctypedef struct base_counts_struct:
    char * base
    int counts

cdef class BaseCounter(object):
    '''
    Class for counting bases in a bam file.
    '''
    cdef Samfile _bam_file
    cdef int min_base_qual
    cdef int min_map_qual

    def __init__(self, Samfile bam_file, int min_base_qual=10, int min_map_qual=10):
        self._bam_file = bam_file
        
        self.min_base_qual = min_base_qual
        self.min_map_qual = min_map_qual
    
    def iter_ref(self, ref):
        '''
        Returns a BaseCounterIter iterator for given ref.
        '''
        iter = BaseCounterIter(ref, self._bam_file.pileup(ref), self.min_base_qual, self.min_map_qual)
        
        return iter
    
    property refs:
        def __get__(self):
            return self._bam_file.references
        
cdef class BaseCounterRow(object):
    '''
    Class for storing count data from Bam file position.
    '''
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
        '''
        1-based position.
        '''    
        def __get__(self):
            return self._position + 1
    
    property depth:
        '''
        Depth of all counts.
        '''
        def __get__(self):
            cdef int depth
            
            depth = self._counts.A + self._counts.C + self._counts.G + self._counts.T
            
            return depth
        
    cdef int get_counts(self, char * base):
        '''
        Lookup the counts for a given base.
        '''
        cdef int counts
    
        counts = 0
        
        if strcmp(base, "A") == 0: 
            counts = self._counts.A
        elif strcmp(base, "C") == 0: 
            counts = self._counts.C
        elif strcmp(base, "G") == 0:
            counts = self._counts.G
        elif strcmp(base, "T") == 0:
            counts = self._counts.T
        
        return counts
    
    cdef base_counts_struct get_base_counts(self, char * base):
        cdef base_counts_struct base_counts
        
        base_counts.base = base
        
        if strcmp(base, "A") == 0: 
            base_counts.counts = self._counts.A
        elif strcmp(base, "C") == 0: 
            base_counts.counts = self._counts.C
        elif strcmp(base, "G") == 0:
            base_counts.counts = self._counts.G
        elif strcmp(base, "T") == 0:
            base_counts.counts = self._counts.T
        
        return base_counts        

cdef class BaseCounterIter(object):
    '''
    Iterator class for iterating over references in bam files and returning BaseCounterRow objects.
    '''
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
        return self.cnext()

    
    cdef BaseCounterRow cnext(self):
        cdef counts_struct counts
        cdef int position
        cdef PileupProxy pileup_column
        cdef BaseCounterRow row
        
        pileup_column = self.pileup_iter.next()
    
        counts = self.parse_pileup_column(pileup_column)
    
        position = pileup_column.pos
    
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
                        
#===============================================================================
# Joint File Analysis
#===============================================================================
cdef class JointBinaryBaseCounter(object):
    '''
    Class for iterating over positions from paired genome files counting bases.
    '''
    cdef BaseCounter _normal_counter
    cdef BaseCounter _tumour_counter
    
    cdef Fastafile _ref_genome_fasta
    
    cdef list _refs
    
    def __init__(self, Samfile normal_bam, Samfile tumour_bam, Fastafile ref_genome_fasta):
        self._normal_counter = BaseCounter(normal_bam)
        self._tumour_counter = BaseCounter(tumour_bam)
        
        self._ref_genome_fasta = ref_genome_fasta
        
        self._refs = list(set(self._normal_counter.refs) & set(self._tumour_counter.refs)) 
        
    def iter_ref(self, ref):
        if ref not in self._refs:
            raise Exception("Invalid reference passed.")
        
        return JointBinaryBaseCounterIterator(
                                               self._normal_counter.iter_ref(ref),
                                               self._tumour_counter.iter_ref(ref),
                                               self._ref_genome_fasta
                                               )        

cdef class JointBinaryCounterRow
cdef JointBinaryCounterRow makeJointBinaryCounterRow(char * ref_base, BaseCounterRow normal_row, BaseCounterRow tumour_row):
    '''
    Constructor method for creating a JointBinaryCounterRow from C.
    '''
    cdef char * non_ref_base

    cdef JointBinaryCounterRow row = JointBinaryCounterRow.__new__(JointBinaryCounterRow)
    
    row._ref = normal_row.ref
    row._position = normal_row.position
    
    row._ref_base = ref_base
     
    # Figure out which is the second most prevalent base in the tumour.
    non_ref_base = get_non_ref_base(tumour_row, ref_base)
    row._non_ref_base = non_ref_base
     
    row._normal_counts = get_binary_counts(ref_base, non_ref_base, normal_row)
    row._tumour_counts = get_binary_counts(ref_base, non_ref_base, tumour_row)    
     
    return row

cdef class JointBinaryCounterRow(object):
    '''
    Class for storing binary count data from a pair of Bam files at a position.
    '''
    cdef char * _ref
    cdef char * _ref_base
    cdef char * _non_ref_base    
    
    cdef int _position
    
    cdef binary_counts_struct _normal_counts
    cdef binary_counts_struct _tumour_counts
    
    def __init__(self):
        raise TypeError("This class cannot be instantiated from Python")
    
    def __str__(self):
        return "\t".join((
                         self._ref,
                         str(self._position),
                         self._ref_base,
                         self._non_ref_base,
                         str(self._normal_counts.A),
                         str(self._normal_counts.B),
                         str(self._tumour_counts.A),
                         str(self._tumour_counts.B)
                         ))
    
    property jcnt_row:
        def __get__(self):
            return [
                    self._position,
                    self._ref_base,
                    self._non_ref_base,
                    self._normal_counts.A,
                    self._normal_counts.B,
                    self._tumour_counts.A,
                    self._tumour_counts.B
                    ]
    
    property ref:
        def __get__(self):
            return self._ref
    
    property position:
        '''
        1-based position.
        '''
        def __get__(self):
            return self._position + 1
    
    property min_depth:
        def __get__(self):
            cdef int normal_depth
            cdef int tumour_depth
            
            normal_depth = self._normal_counts.A + self._normal_counts.B
            tumour_depth = self._tumour_counts.A + self._tumour_counts.B
            
            if normal_depth < tumour_depth:
                return normal_depth
            else:
                return tumour_depth
    
    property has_var:
        def __get__(self):
            if self._tumour_counts.B > 0:
                return True
            else:
                return False
            
cdef class JointBinaryBaseCounterIterator(object):
    cdef BaseCounterIter _normal_iter
    cdef BaseCounterIter _tumour_iter
    
    cdef Fastafile _ref_genome_fasta
    
    cdef BaseCounterRow _normal_row
    cdef BaseCounterRow _tumour_row

    def __init__(self, BaseCounterIter normal_iter, BaseCounterIter tumour_iter, Fastafile ref_genome_fasta):
        self._normal_iter = normal_iter
        self._tumour_iter = tumour_iter
        
        self._ref_genome_fasta = ref_genome_fasta
        
        try:
            self._normal_row = normal_iter.cnext()
            self._tumour_row = tumour_iter.cnext()
        except StopIteration:
            raise Exception("Empty iterator passed.")
    
    def __iter__(self):
        return self
    
    def __next__(self):
        cdef JointBinaryCounterRow bin_count_row
        
        cdef BaseCounterRow normal_row
        cdef BaseCounterRow tumour_row
        
        cdef int normal_pos
        cdef int tumour_pos
        
        while True:
            normal_pos = self._normal_row._position
            tumour_pos = self._tumour_row._position
            
            if normal_pos == tumour_pos:
                bin_count_row = self._get_joint_binary_counter_row()
                
                self._normal_row = self._normal_iter.cnext()
                self._tumour_row = self._tumour_iter.cnext()
                
                return bin_count_row
            elif normal_pos < tumour_pos:
                self._normal_row = self._normal_iter.cnext()
            elif normal_pos > tumour_pos:
                self._tumour_row = self._tumour_iter.cnext()
            else:
                raise Exception("Error in joint pileup iterator.")
    
    cdef JointBinaryCounterRow _get_joint_binary_counter_row(self):
        cdef int region_length
        cdef char * ref_base
        cdef JointBinaryCounterRow bin_count_row
        
        region_length = 1 
        
        ref_base = self._ref_genome_fasta._fetch(
                                                 self._normal_row._ref,
                                                 self._normal_row._position,
                                                 self._normal_row._position + 1,
                                                 & region_length
                                                 )
    
        bin_count_row = makeJointBinaryCounterRow(ref_base, self._normal_row, self._tumour_row)
        
        return bin_count_row            
                        
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

#===============================================================================
# Utility functions for finding non-ref bases and counts.
#===============================================================================
cdef int compare_base_counts_struct(const_void * a, const_void * b):
    cdef base_counts_struct * first
    cdef base_counts_struct * second
    
    first = < base_counts_struct *> a
    second = < base_counts_struct *> b

    return first.counts - second.counts

cdef char * get_non_ref_base(BaseCounterRow row, char * ref_base):
    '''
    Sort the bases by number observed and return the most common non-reference base.
    ''' 
    cdef base_counts_struct[4] counts

    counts[0] = row.get_base_counts("A")
    counts[1] = row.get_base_counts("C")
    counts[2] = row.get_base_counts("G")
    counts[3] = row.get_base_counts("T")
    
    qsort(counts, 4, sizeof(int), compare_base_counts_struct)
    
    # If the most prevalent base is not the reference return it.
    if strcmp(counts[3].base, ref_base) != 0:
        return counts[3].base
    # Otherwise return the second most prevalent base.
    else:
        return counts[2].base

cdef binary_counts_struct get_binary_counts(char * ref_base, char * non_ref_base, BaseCounterRow row):
    cdef binary_counts_struct counts

    counts.A = row.get_counts(ref_base)
    counts.B = row.get_counts(non_ref_base)
    
    return counts 
