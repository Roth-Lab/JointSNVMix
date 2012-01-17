'''
Created on 2012-01-17

@author: Andrew Roth
'''
cdef class BamFile:
    cdef bam_index_t * _index

    def __cinit__(self, char * file_name):
        self._file_name = file_name
        
        self._open_file()
        self._load_index_file()
    
    def __dealloc__(self):
        bam_index_destroy(self._index)

    cdef _load_index_file(self, file_name):
        index_file_name = self._file_name + ".bai"
        
        if not os.path.exists(index_file_name):
            raise Exception("Index not found for BAM file {0}".format(file_name))
        
        self._index = bam_index_load(file_name)
    
    cdef _open_file(self):
        self._bam_file = samopen(filename, mode, NULL)
        
        if self._bam_file == NULL:
            raise Exception("Could not open file {0} - is it SAM/BAM format?".format(self._file_name))

        if self._bam_file.header == NULL:
            raise Exception("File {0} does not have valid header - is it SAM/BAM format?".format(self._file_name))
    
    cdef  samfile_t * get_file_pointer(self):
        return self._bam_file
        
    def pileup(self, reference=None, start=None, end=None, **kwargs):
        cdef int rtid, rstart, rend, has_coord
        cdef bam_plbuf_t * buf

        has_coord, rtid, rstart, rend = self._parseRegion(reference, start, end, region)


        if has_coord:
            return IteratorColumnRegion(self,
                                         tid=rtid,
                                         start=rstart,
                                         end=rend,
                                         **kwargs)
        else:
            return IteratorColumnAllRefs(self, **kwargs)
        
    def _parseRegion(self, reference=None, start=None, end=None, region=None):
        '''
        Parse region information.

        raise ValueError for for invalid regions.

        returns a tuple of flag, tid, start and end. Flag indicates
        whether some coordinates were supplied.

        Note that regions are 1-based, while start,end are python coordinates.
        '''
        # This method's main objective is to translate from a reference to a tid. 
        # For now, it calls bam_parse_region, which is clumsy. Might be worth
        # implementing it all in pysam (makes use of khash).
        
        cdef int rtid
        cdef long long rstart
        cdef long long rend

        rtid = -1
        rstart = 0
        rend = max_pos
        
        if start != None: 
            try:
                rstart = start
            except OverflowError:
                raise ValueError('start out of range (%i)' % start)
            
        if end != None: 
            try:
                rend = end
            except OverflowError:
                raise ValueError('end out of range (%i)' % end)

        if region:
            parts = re.split("[:-]", region)
            reference = parts[0]
            if len(parts) >= 2: rstart = int(parts[1]) - 1
            if len(parts) >= 3: rend = int(parts[2])

        if not reference: return 0, 0, 0, 0

        rtid = self.gettid(reference)
        if rtid < 0: raise ValueError("invalid reference `%s`" % reference)
        if rstart > rend: raise ValueError('invalid coordinates: start (%i) > end (%i)' % (rstart, rend))
        if not 0 <= rstart < max_pos: raise ValueError('start out of range (%i)' % rstart)
        if not 0 <= rend <= max_pos: raise ValueError('end out of range (%i)' % rend)

        return 1, rtid, rstart, rend

cdef class AlignedRead:
    def __dealloc__(self):
        bam_destroy1(self._delegate)

    def overlap(self, uint32_t start, uint32_t end):
        """return number of aligned bases of read overlapping the interval *start* and *end*
        on the reference sequence.
        """
        cdef uint32_t k, i, pos, overlap
        cdef int op, o
        cdef uint32_t * cigar_p
        cdef bam1_t * src 

        overlap = 0

        src = self._delegate
        if src.core.n_cigar == 0: return 0
        pos = src.core.pos
        o = 0

        cigar_p = bam1_cigar(src)
        for k from 0 <= k < src.core.n_cigar:
            op = cigar_p[k] & BAM_CIGAR_MASK
            l = cigar_p[k] >> BAM_CIGAR_SHIFT

            if op == BAM_CMATCH:
                o = min(pos + l, end) - max(pos, start)
                if o > 0: overlap += o

            if op == BAM_CMATCH or op == BAM_CDEL or op == BAM_CREF_SKIP:
                pos += l

        return overlap

    def opt(self, tag):
        """retrieves optional data given a two-letter *tag*"""
        #see bam_aux.c: bam_aux_get() and bam_aux2i() etc 
        cdef uint8_t * v
        v = bam_aux_get(self._delegate, tag)
        if v == NULL: raise KeyError("tag '%s' not present" % tag)
        type = chr(v[0])
        if type == 'c' or type == 'C' or type == 's' or type == 'S':
            return < int > bam_aux2i(v)            
        elif type == 'i' or type == 'I':
            return < int32_t > bam_aux2i(v)            
        elif type == 'f' or type == 'F':
            return < float > bam_aux2f(v)
        elif type == 'd' or type == 'D':
            return < double > bam_aux2d(v)
        elif type == 'A':
            # there might a more efficient way
            # to convert a char into a string
            return '%c' % < char > bam_aux2A(v)
        elif type == 'Z':
            return < char *> bam_aux2Z(v)


cdef makeAlignedRead(bam1_t * src):
    cdef AlignedRead dest = AlignedRead.__new__(AlignedRead)
    
    dest._delegate = bam_dup1(src)
    
    return dest

cdef class IteratorRowRegion(object):
    cdef int _return_value    
    
    cdef bam_iter_t _iter
    cdef bam1_t * _bam_struct
    cdef samfile_t * _bam_file_ptr
    
    cdef BamFile _bam_file

    def __cinit__(self, BamFile bam_file, int tid, int start, int end, int reopen=True):        
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

    cdef int cnext(self):
        self._return_value = bam_iter_read(self._bam_file_ptr.x.bam,
                                           self._iter,
                                           self._bam_struct)
        
    def __next__(self): 
        self.cnext()

        if self.retval < 0: 
            raise StopIteration
        
        return makeAlignedRead(self._bam_struct)

cdef class IterData(object):
    cdef char * _seq
    
    cdef int _seq_len
    cdef int _tid
    
    cdef samfile_t * _bam_file_ptr
    
    cdef faidx_t * _fasta_file_ptr
    
    cdef bam_iter_t _iter

    def __cinit__(self, samfile_t * bam_file_ptr, faidx_t * fasta_file_ptr, bam_iter_t iter):        
        self._bam_file_ptr = bam_file_ptr
        
        self._fasta_file_ptr = fasta_file_ptr
        
        self._iter = iter
        
        self._seq = NULL
        self._tid = -1

cdef class IteratorColumn:
    cdef int _tid
    cdef int _pos
    cdef int _n_plp
    cdef int _mask
    
    cdef BamFile _bam_file
    cdef FastaFile _fasta_file
    
    cdef const_bam_pileup1_t_ptr _plp    
    cdef bam_plp_t _pileup_iter    
    cdef __iterdata _iter_data 

    def __cinit__(self, BamFile bam_file, FastaFile fasta_file):        
        self._bam_file = bam_file
        
        self._fasta_file = fasta_file
                                
        self._mask = BAM_DEF_MASK

        self._tid = 0
        self._pos = 0
        self._n_plp = 0
        self._plp = NULL
                
        self._setup_iterator_data(tid, start, end, 1)

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
            
            return makePileupProxy(< bam_pileup1_t *> self.plp,
                                     self.tid,
                                     self.pos,
                                     self.n_plp)    
    
    cdef int cnext(self):
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
        
        self._iter_data.bam_file_prt = self._bam_file.get_file_ptr()
        self._iter_data.fasta_file_prt = self._fasta_file.get_file_ptr()
        self._iter_data.iter = self._iter
        self._iter_data.seq = NULL
        self._iter_data.tid = -1
            
        self._pileup_iter = bam_plp_init(& __advance_all, & self._iter_data)
    
        bam_plp_set_mask(self._pileup_iter, self._mask)
        
#=======================================================================================================================
# Utility functions and structs
#=======================================================================================================================
ctypedef struct __iter_data:
    int tid
    int seq_len
    
    char * seq
    
    bam_iter_t iter
    
    faidx_t * fasta_file_ptr
    
    samfile_t * bam_file_ptr

cdef int __advance_all(void * data, bam1_t * b):
    '''
    Advance without any read filtering.
    '''
    cdef __iter_data * d
    
    d = < __iterdata *> data
    
    return bam_iter_read(d.bam_file_ptr.x.bam, d.iter, b)
