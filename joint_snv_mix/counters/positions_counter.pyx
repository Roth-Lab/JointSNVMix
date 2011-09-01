'''
Created on 2011-08-11

@author: Andrew Roth
'''
cdef class PositionsCounter(Counter):
    def __init__(self, Counter counter, char * pos_file_name):
        self._starts = []
        self._pos_file_name = pos_file_name
        self._counter = counter
        
        self._load_index()
        
        # Restrict to refs present if both positions file and counter
        self._refs = tuple(set(self.refs) & set(counter.refs))
    
    property refs:
        def __get__(self):
            return self._refs
    
    def iter_ref(self, ref):
        if ref not in self.refs:
            raise Exception('Invalid reference.')

        i = self._refs.index(ref)
        start = self._starts[i]
        stop = self._starts[i + 1]
        
        pos_iter = PositionsIterator(self._pos_file_name, start, stop)
        
        counter_iter = self._counter.iter_ref(ref)
        
        return PositionsCounterRefIterator(ref, counter_iter, pos_iter)
    
    cdef _load_index(self):
        cdef FILE * file_p
        cdef char ref[100], prev_ref[100]
        cdef int pos, coord
        cdef int ref_result, coord_result
        cdef list refs = []
        
        prev_ref[0] = '\0'
        
        file_p = fopen(self._pos_file_name, "r")
        
        if file_p == NULL:
            raise Exception("Couldn't open positions file")
        
        pos = 0
        
        while True:
            ref_result = fscanf(file_p, "%s" , ref)
            coord_result = fscanf(file_p, "%d", & coord)
            
            if ref_result == EOF:
                break                        
            
            if strcmp(prev_ref, ref) != 0:
                refs.append(ref)
                self._starts.append(pos)
                
                strcpy(prev_ref, ref)
            
            pos = ftell(file_p)
        
        self._starts.append(pos)
        
        fclose(file_p)
        
        self._refs = tuple(refs)
        
cdef class PositionsCounterRefIterator(RefIterator):
    def __init__(self, char * ref, RefIterator ref_iter, PositionsIterator pos_iter):
        self._ref = ref
        self._ref_iter = ref_iter
        self._pos_iter = pos_iter
        
    cdef advance_position(self):
        cdef int pos_1, pos_2
        
        self._ref_iter.advance_position()
        self._pos_iter.cnext()
        
        while True:
            pos_1 = self._ref_iter._position
            pos_2 = self._pos_iter._position

            if pos_1 == pos_2:
                self._position = pos_1
                
                break
            elif pos_1 < pos_2:
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
            raise Exception("Couldn't open the spam file")
        
        fseek(self._file_p, start, 0)
    
    def __iter__(self):
        return self
    
    def __next__(self):
        self.cnext()
        
        return self.position
    
    property position:
        def __get__(self):
            return self._position + 1
    
    cdef cnext(self):
        cdef char ref[100]
        cdef int coord
        cdef int ref_result, coord_result
        
        if ftell(self._file_p) >= self._stop:            
            raise StopIteration
        
        ref_result = fscanf(self._file_p, "%s" , ref)
        coord_result = fscanf(self._file_p, "%d", & coord)
        
        # Convert to 0-based
        self._position = coord - 1
    
    def __dealloc__(self):
        fclose(self._file_p)
