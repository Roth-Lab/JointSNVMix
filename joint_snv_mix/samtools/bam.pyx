'''
Created on 2012-01-17

@author: Andrew Roth
'''
import os

# Global constant
cdef int max_pos = 2 << 29

cdef class BamFile:
    '''
    Heavily simplified version of pysam Samfile class. Tailored specifically for extracting data from BAM files for
    JointSNVMix.
    '''
    def __cinit__(self, char * file_name):
        self._file_name = strdup(file_name)
        
        self._open_file()

        self._load_index()
        
        self._load_references()
            
    def __dealloc__(self):
        if self._file_name != NULL:
            free(self._file_name)
        
        if self._index != NULL:
            bam_index_destroy(self._index)
        
        if self._bam_file != NULL:
            samclose(self._bam_file)
    
    def _open_file(self):
        self._bam_file = samopen(self._file_name, 'rb', NULL)
        
        if self._bam_file == NULL:
            raise Exception("Could not open file {0} - is it SAM/BAM format?".format(self._file_name))

        if self._bam_file.header == NULL:
            raise Exception("File {0} does not have valid header - is it SAM/BAM format?".format(self._file_name))        
        
    def _load_index(self):
        self._index = bam_index_load(self._file_name)
        
        if self._index == NULL:
            raise Exception("Index not found for BAM file {0}".format(self._file_name))
    
    def _load_references(self):        
        self._refs = []
        self._tids = {}
            
        for i in range(self._bam_file.header.n_targets):
            self._refs.append(self._bam_file.header.target_name[i])
        
        for ref in self._refs:
            tid = self._get_tid(ref)

            self._tids[ref] = tid
            
    def get_pileup_iterator(self, reference, start=None, stop=None):
        if reference not in self._refs:
            raise Exception("Invalid reference given `{0}`".format(reference))
        
        tid = self._tids[reference]

        if start is not None and stop is not None:
            return PileupIterator(self, tid=tid, start=start, stop=stop)
        else:
            return PileupIterator(self, tid=tid) 

    def _get_tid(self, reference):
        '''
        Use bam_parse_region to extract the tid, start and end of a reference.
        '''
        cdef int tid, start, end

        bam_parse_region(self._bam_file.header, reference, & tid, & start, & end)        
        
        if tid < 0: 
            raise Exception('BamFile.pileup : Invalid reference {0}.'.format(reference))

        return tid    
    
    cdef samfile_t * get_file_pointer(self):
        return self._bam_file
    
    cdef bam_index_t * get_index(self):
        return self._index

    property references:
        def __get__(self):
            return tuple(self._refs)
