'''
Created on 2012-01-25

@author: Andrew Roth
'''
cdef class PositionsCounter(object):
    '''
    Class for iterating only positions specified in file.
    '''
    def __init__(self, char * pos_file_name, counter):
        self._pos_file_name = pos_file_name
        
        self._load_index()
        
        self._refs = tuple(self._index.keys())
        
        # Restrict to refs present if both positions file and counter
        self._refs = tuple(set(self.refs) & set(counter.refs))
    
    def get_ref_iterator(self, counter_iter):
        ref = counter_iter.ref
        
        if ref not in self.refs:
            raise Exception('Invalid reference.')

        start, stop = self._index[ref]
        
        pos_iter = PositionsIterator(self._pos_file_name, start, stop)
        
        return PositionsCounterRefIterator(ref, counter_iter, pos_iter)
    
    cdef _load_index(self):
        '''
        For efficiency index the first and last lines of a file where each chromosome appears.
        '''
        cdef FILE * file_p
        cdef char ref[100], prev_ref[100]
        cdef int file_pos, ref_start, coord
        
        self._index = {}
        
        prev_ref[0] = '\0'
        
        file_p = fopen(self._pos_file_name, "r")
        
        if file_p == NULL:
            raise Exception("Couldn't open the positions file {0}".format(self._pos_file_name))
        
        file_pos = 0
        ref_start = file_pos
        
        while True:
            if fscanf(file_p, "%s" , ref) == EOF:
                break
            
            fscanf(file_p, "%d", & coord)                    
            
            if strcmp(prev_ref, ref) != 0:
                if file_pos > 0:
                    self._index[prev_ref] = (ref_start, file_pos)
                    ref_start = file_pos
                                              
                strcpy(prev_ref, ref)
            
            file_pos = ftell(file_p)
        
        self._index[ref] = (ref_start, file_pos)
        
        fclose(file_p)        
        
    property refs:
        def __get__(self):
            return self._refs        
        
cdef class PositionsCounterRefIterator(RefIterator):
    def __init__(self, char * ref, JointBinaryCounterIterator ref_iter, PositionsIterator pos_iter):
        self._ref = ref
        self._ref_iter = ref_iter
        self._pos_iter = pos_iter
        
    cdef advance_position(self):
        cdef int pos_1, pos_2
        
        self._ref_iter.advance_position()
        self._pos_iter.cnext()
        
        while True:
            pos_1 = self._ref_iter._pos
            pos_2 = self._pos_iter._pos

            if pos_1 == pos_2:
                self._pos = pos_1
                
                break
            elif pos_1 < pos_2:
                # This could cause problems it the length of reads exceed 100000.
                if pos_2 - pos_1 > 100000:
                    self._ref_iter.jump_to_position(pos_2)
                else:                
                    self._ref_iter.advance_position()                        
            elif pos_1 > pos_2:
                self._pos_iter.cnext()                
            else:
                raise Exception("Invalid iteration occured.")
    
    cdef parse_current_position(self):
        self._ref_iter.parse_current_position()
        
        self._current_row = self._ref_iter._current_row                 

cdef class PositionsIterator(object):
    def __init__(self, char * pos_file_name, int start, int stop):
        self._stop = stop
        
        self._file_p = fopen(pos_file_name, "r")
        
        if self._file_p == NULL:
            raise Exception("Couldn't open the positions file {0}".format(pos_file_name))
        
        fseek(self._file_p, start, 0)

    def __dealloc__(self):
        fclose(self._file_p)
    
    cdef cnext(self):
        cdef char ref[100]
        cdef int coord
        cdef int ref_result, coord_result
        
        if ftell(self._file_p) >= self._stop:            
            raise StopIteration
        
        ref_result = fscanf(self._file_p, "%s" , ref)
        coord_result = fscanf(self._file_p, "%d", & coord)
        
        # Convert to 0-based
        self._pos = coord - 1
