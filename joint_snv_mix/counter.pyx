'''
Classes for creating iterators for count data over a pair of genomes. 

Created on 2012-01-18

@author: Andrew Roth
'''
cdef class JointBinaryCounter(object):
    '''
    Class for iterating over positions from paired genome files counting bases.
    
    Parameters :
        normal_counter - A BAM file object for the normal genome.
        tumour_counter - A BAM file object for the tumour genome.
        ref_genome - A FastaFile for the reference genome.
    '''
    def __init__(self,
                 BamFile normal_bam,
                 BamFile tumour_bam,
                 FastaFile ref_genome,
                 int min_base_qual=0,
                 int min_map_qual=0,
                 bint qualities=0):
        
        self._qualities = qualities
        
        self._min_base_qual = min_base_qual
        self._min_map_qual = min_map_qual
        
        self._normal_bam = normal_bam
        self._tumour_bam = tumour_bam
        
        self._ref_genome = ref_genome
        
        self._refs = tuple(set(self._normal_bam.references) & set(self._tumour_bam.references))
    
    def get_ref_iterator(self, ref):
        if ref not in self.refs:
            raise Exception("Invalid reference passed.")
        
        return JointBinaryCounterIterator(ref,
                                          self._normal_bam.get_pileup_iterator(ref),
                                          self._tumour_bam.get_pileup_iterator(ref),
                                          self._ref_genome,
                                          self._min_base_qual,
                                          self._min_map_qual,
                                          self._qualities)

    property refs:
        '''
        Read only access to list of available references.
        '''
        def __get__(self):
            return self._refs        
        
cdef class JointBinaryCounterIterator(RefIterator):
    def __init__(self,
                 char * ref,
                 PileupIterator normal_iter,
                 PileupIterator tumour_iter,
                 FastaFile ref_genome,
                 int min_base_qual,
                 int min_map_qual,
                 bint qualities):
              
        self._ref = ref
        self._pos = -1
        
        self._normal_iter = normal_iter
        self._tumour_iter = tumour_iter

        self._row_factory = RowFactory(ref_genome, min_base_qual, min_map_qual, qualities)

    cdef advance_position(self):        
        cdef int normal_pos
        cdef int tumour_pos
        
        self._normal_iter.advance_position()        
        self._tumour_iter.advance_position()
                
        while True:
            normal_pos = self._normal_iter._pos
            tumour_pos = self._tumour_iter._pos
            
            if normal_pos == tumour_pos:
                self._pos = normal_pos                                       
                break
            elif normal_pos < tumour_pos:
                self._normal_iter.advance_position()
            elif normal_pos > tumour_pos:
                self._tumour_iter.advance_position()
            else:
                raise Exception("Error in joint pileup iterator.")        
    
    cdef parse_current_position(self):
        cdef char * ref_base
        cdef PileupColumn normal_column
        cdef PileupColumn tumour_column
        
        self._normal_iter.parse_current_position()        
        normal_column = self._normal_iter._column
        
        self._tumour_iter.parse_current_position()
        tumour_column = self._tumour_iter._column
                
        self._current_row = self._row_factory.get_row(self._ref, self._pos, normal_column, tumour_column)
    
    cdef jump_to_position(self, int position):
        cdef int normal_pos
        cdef int tumour_pos
    
        self._normal_iter.jump_to_position(position)
        self._tumour_iter.jump_to_position(position)
        
        
        self.advance_position()

        normal_pos = self._normal_iter._pos
        tumour_pos = self._tumour_iter._pos
        
        if normal_pos != tumour_pos:
            raise Exception('Jump did not work.')
        else:
            self._pos = normal_pos
