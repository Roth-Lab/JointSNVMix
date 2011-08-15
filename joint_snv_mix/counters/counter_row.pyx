'''
Created on 2011-07-29

@author: Andrew Roth
'''
cdef class CounterRow(object):
    def __init__(self):
        raise TypeError("This class cannot be instantiated from Python")
    
    property ref:
        def __get__(self):
            return self._ref
    
    property position:
        '''
        1-based position
        '''
        def __get__(self):
            return self._position + 1

cdef class SingleSampleCounterRow(CounterRow):    
    property depth:
        def __get__(self):    
            return self._depth
    
    cdef int get_counts(self, char * base):
        '''
        Return the counts for given base.
        '''
        pass

cdef class PairedSampleBinomialCounterRow(CounterRow):
    def __str__(self):
        out_row = [self.ref, str(self.position), self.ref_base, self.non_ref_base]
        out_row.extend([str(x) for x in self.counts])
        
        return "\t".join(out_row)
    
    property counts:
        def __get__(self):
            return [x for x in self._counts[:4]]

    property normal_depth:
        def __get__(self):
            return self._normal_depth
    
    property tumour_depth:
        def __get__(self):
            return self._tumour_depth
    
    property ref_base:
        def __get__(self):
            return self._ref_base
        
    property non_ref_base:
        def __get__(self):
            return self._non_ref_base
#===============================================================================
# Utility functions for finding non-ref bases and counts.
#===============================================================================
cdef int compare_base_counts_struct(const_void * a, const_void * b):
    cdef base_counts_struct * first
    cdef base_counts_struct * second
    
    first = < base_counts_struct *> a
    second = < base_counts_struct *> b

    return first.counts - second.counts

cdef char * get_non_ref_base(char * ref_base, SingleSampleCounterRow normal_row, SingleSampleCounterRow tumour_row):
    '''
    Sort the bases by number observed and return the most common non-reference base.
    
    Return N if the counts of most common non-reference is 0.
    ''' 
    cdef int non_ref_index    
    cdef base_counts_struct[4] counts

    counts[0].base = "A"
    counts[1].base = "C"
    counts[2].base = "G"
    counts[3].base = "T"
    
    counts[0].counts = normal_row.get_counts("A") + tumour_row.get_counts("A")
    counts[1].counts = normal_row.get_counts("C") + tumour_row.get_counts("C")
    counts[2].counts = normal_row.get_counts("G") + tumour_row.get_counts("G")
    counts[3].counts = normal_row.get_counts("T") + tumour_row.get_counts("T")

    # Sort the structs by counts field. 
    qsort(counts, 4, sizeof(base_counts_struct), compare_base_counts_struct)
            
    # If the most prevalent base is not the reference return it.
    if strcmp(counts[3].base, ref_base) != 0:
        non_ref_index = 3
    else:
        non_ref_index = 2
    
    if counts[non_ref_index].counts > 0:
        return counts[non_ref_index].base
    else:
        return 'N'
        
cdef binary_counts_struct get_binary_counts(char * ref_base, char * non_ref_base, SingleSampleCounterRow row):
    cdef binary_counts_struct counts

    counts.A = row.get_counts(ref_base)
    counts.B = row.get_counts(non_ref_base)
    
    return counts        
