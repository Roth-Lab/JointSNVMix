'''
Created on 2012-01-17

@author: Andrew Roth
'''

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