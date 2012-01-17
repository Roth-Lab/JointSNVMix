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

        self._load_index(file_name)
            
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
            raise Exception("Index not found for BAM file {0}".format(file_name))
        
        

    def pileup(self, reference, start=None, end=None):
        cdef int region_tid, region_start, region_end

        region_tid, region_start, region_end = self._parseRegion(reference, start, end)

        return IteratorColumnRegion(self, tid=region_tid, start=region_start, end=region_end)    

    def _parseRegion(self, reference, start=None, end=None):
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
    
    cdef samfile_t * get_file_pointer(self):
        return self._bam_file
    
    cdef bam_index_t * get_index(self):
        return self._index

    property references:
        def __get__(self):
            t = []
            
            for x in range(self._bam_file.header.n_targets):
                t.append(self._bam_file.header.target_name[x])
            
            return tuple(t)
