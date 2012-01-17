'''
A simple wrapper around a core set of samtools functionality. Based on code from the pysam library. Heavily modified, to
simplify access for JointSNVMix.

Created on 2012-01-17

@author: Andrew Roth
'''
import os

#=======================================================================================================================
# Global constants
#=======================================================================================================================
cdef int max_pos = 2 << 29

#=======================================================================================================================
# Fasta file
#=======================================================================================================================
cdef class FastaFile:
    def __cinit__(self, char * file_name):
        self._file_name = strdup(file_name)
        
        self._fasta_file = fai_load(file_name)

        if self._fasta_file == NULL:
            raise Exception("Could not open FASTA file {0}".format(file_name))
        
        index_file_name = file_name + "fai"
        
        if not os.path.exists(index_file_name):
            raise Exception("Index not found for FASTA file {0}".format(file_name))

    def close(self):
        if self._fasta_file != NULL:
            fai_destroy(self._fasta_file)            
            self._fasta_file = NULL

    def __dealloc__(self):
        self.close()
        
        if self._file_name != NULL: 
            free(self._file_name)
    
    cdef faidx_t * get_file_pointer(self):
        return self._fasta_file

    cdef char * get_reference_base(self, char * reference, int position):
        cdef char * ref_base
        cdef int len
        
        len = 1
        
        ref_base = self._fetch(reference, position, position + 1, & len)
        
        ref_base[0] = < char > toupper(< int > ref_base[0])

        return ref_base    

    cdef char * _fetch(self, char * reference, int start, int end, int * length):
        return faidx_fetch_seq(self._fasta_file, reference, start, end - 1, length)

#=======================================================================================================================
# BAM file class
#=======================================================================================================================
cdef class BamFile:
    '''
    Heavily simplified version of pysam Samfile class. Tailored specifically for extracting data from BAM files for
    JointSNVMix.
    '''
    def __cinit__(self, char * file_name):
        self._file_name = file_name
        
        self._bam_file = samopen(file_name, 'rb', NULL)
        
        if self._bam_file == NULL:
            raise Exception("Could not open file {0} - is it SAM/BAM format?".format(self._file_name))

        if self._bam_file.header == NULL:
            raise Exception("File {0} does not have valid header - is it SAM/BAM format?".format(self._file_name))

        self._load_index(file_name)
            
    def __dealloc__(self):
        bam_index_destroy(self._index)
        
    def _load_index(self, file_name):
        index_file_name = self._file_name + ".bai"
        
        if not os.path.exists(index_file_name):
            raise Exception("Index not found for BAM file {0}".format(file_name))
        
        self._index = bam_index_load(file_name)
    
    cdef samfile_t * get_file_pointer(self):
        return self._bam_file
    
    cdef bam_index_t * get_index(self):
        return self._index
        
    def pileup(self, reference, start=None, end=None):
        cdef int region_tid, region_start, region_end

        region_tid, region_start, region_end = self._parseRegion(reference, start, end)

        return IteratorColumnRegion(self, tid=region_tid, start=region_start, end=region_end)

    def _parseRegion(self, reference=None, start=None, end=None):
        '''
        Parse region information.

        raise ValueError for for invalid regions.

        returns a tuple of flag, tid, start and end. Flag indicates
        whether some coordinates were supplied.

        Note that regions are 1-based, while start,end are python coordinates.
        '''
        cdef int rtid
        cdef int rstart
        cdef int rend

        rtid = -1
        rstart = 0
        rend = max_pos
        
        if start != None: 
            try:
                rstart = start
            except OverflowError:
                raise Exception('BamFile.pileup : Start out of numerical range {0}'.format(start))
            
        if end != None: 
            try:
                rend = end
            except OverflowError:
                raise Exception('BamFile.pileup : End out of numerical range {0}'.format(end))

        if start != None and end != None:
            region = "{0}:{1}-{2}".format(reference, start + 1, end)
        else:
            region = reference

        bam_parse_region(self._bam_file.header, region, & rtid, & rstart, & rend)        
        
        if rtid < 0: 
            raise Exception('BamFile.pileup : Invalid reference {0}.'.format(reference))
        
        if rstart > rend: 
            raise Exception('BamFile.pileup : Invalid coordinates: start {0} > end {1}.'.format(rstart, rend))
        
        if not 0 <= rstart < max_pos: 
            raise Exception('BamFile.pileup : Start out of range {0}'.format(rstart))
        
        if not 0 <= rend <= max_pos: 
            raise ValueError('BamFile.pileup : End out of range {0}'.format(rend))

        return rtid, rstart, rend

    property references:
        def __get__(self):
            t = []
            
            for x in range(self._bam_file.header.n_targets):
                t.append(self._bam_file.header.target_name[x])
            
            return tuple(t)    

#=======================================================================================================================
# Aligned Read
#=======================================================================================================================
cdef class AlignedRead:
    def __dealloc__(self):
        bam_destroy1(self._delegate)

cdef makeAlignedRead(bam1_t * src):
    cdef AlignedRead dest = AlignedRead.__new__(AlignedRead)
    
    dest._delegate = bam_dup1(src)
    
    return dest

#=======================================================================================================================
# Iterators
#=======================================================================================================================
cdef class IteratorRowRegion(object):
    def __cinit__(self, BamFile bam_file, int tid, int start, int end):        
        # Makes sure that samfile stays alive as long as the iterator is alive.
        self._bam_file = bam_file

        self._bam_file_ptr = self._bam_file.get_file_pointer()

        self._return_value = 0

        self._iter = bam_iter_query(self._bam_file.get_index(), tid, start, end)

        self._bam_struct = bam_init1()
        
    def __dealloc__(self):
        bam_destroy1(self._bam_struct)        

    def __iter__(self):
        return self 

    cdef bam1_t * getCurrent(self):
        return self._bam_struct

    cdef cnext(self):
        self._return_value = bam_iter_read(self._bam_file_ptr.x.bam,
                                           self._iter,
                                           self._bam_struct)
        
    def __next__(self): 
        self.cnext()

        if self.retval < 0: 
            raise StopIteration
        
        return makeAlignedRead(self._bam_struct)

cdef class IteratorColumnRegion:
    def __cinit__(self, BamFile bam_file, int tid=0, int start=0, int end=max_pos):        
        self._bam_file = bam_file
        
#        self._fasta_file = fasta_file
                                
        self._mask = BAM_DEF_MASK

        self._tid = 0
        self._pos = 0
        self._n_plp = 0
        self._plp = NULL
                
        self._setup_iterator_data(self._tid, start, end)

    def __dealloc__(self):
        # reset in order to avoid memory leak messages for iterators that have
        # not been fully consumed
        if self._pileup_iter != < bam_plp_t > NULL:
            bam_plp_reset(self._pileup_iter)
            
            bam_plp_destroy(self._pileup_iter)
            
            self._pileup_iter = < bam_plp_t > NULL
    
        if self._iter_data.seq != NULL: 
            free(self._iter_data.seq)
            
            self._iter_data.seq = NULL
        
    def __iter__(self):
        return self

    def __next__(self): 
        while True:
            self.cnext()
            
            if self._n_plp < 0:
                raise Exception("IteratorColumn : error during iteration.")
        
            if self._plp == NULL:
                raise StopIteration
            
            return makePileupProxy(< bam_pileup1_t *> self._plp,
                                     self._tid,
                                     self._pos,
                                     self._n_plp)    
    
    cdef cnext(self):
        '''
        Perform next iteration.
    
        This method is analogous to the samtools bam_plp_auto method.
        It has been re-implemented to permit for filtering.
        '''
        self._plp = bam_plp_auto(self._pileup_iter,
                                 & self._tid,
                                 & self._pos,
                                 & self._n_plp)
        
    cdef _setup_iterator_data(self, int tid, int start, int end):
        '''
        Setup the iterator structure.
        '''    
        self._iter = IteratorRowRegion(self._bam_file, tid, start, end)
        
        self._iter_data.bam_file_ptr = self._bam_file.get_file_pointer()
        self._iter_data.fasta_file_ptr = NULL
        self._iter_data.iter = self._iter._iter
        self._iter_data.seq = NULL
        self._iter_data.tid = -1
            
        self._pileup_iter = bam_plp_init(& __advance_all, & self._iter_data)
    
        bam_plp_set_mask(self._pileup_iter, self._mask)

#=======================================================================================================================
# Pileup objects
#=======================================================================================================================
cdef class PileupProxy:
    pass

cdef class PileupRead:
    pass

# Factory methods
#---------------------------------------------------------------------------------------------------------------------- 
cdef makePileupProxy(bam_pileup1_t * plp, int tid, int pos, int n):
    cdef PileupProxy dest = PileupProxy.__new__(PileupProxy)
    
    dest._plp = plp
    dest._tid = tid
    dest._pos = pos
    dest._n_pu = n
    
    return dest

cdef makePileupRead(bam_pileup1_t * src):
    cdef PileupRead dest = PileupRead.__new__(PileupRead)
    
    dest._alignment = makeAlignedRead(src.b)
    dest._qpos = src.qpos
    dest._indel = src.indel
    dest._level = src.level
    dest._is_del = src.is_del
    dest._is_head = src.is_head
    dest._is_tail = src.is_tail
    
    return dest
        
#=======================================================================================================================
# Utility functions and structs
#=======================================================================================================================
cdef int __advance_all(void * data, bam1_t * b):
    '''
    Advance without any read filtering.
    '''
    cdef __iter_data * d
    
    d = < __iter_data *> data
    
    return bam_iter_read(d.bam_file_ptr.x.bam, d.iter, b)
