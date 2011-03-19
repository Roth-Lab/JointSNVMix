'''
Abstract bases classes that all file format IO objects should inherit from. The goal is to provide a consistent
interface for reading and writing file objects stored in HDF5 format. The underlying file structures should be hidden
completely from the interface.
'''
class BinaryReader(object):
    '''
    Class for reading BinaryFile object.
    
    Should be able to iterate over the data group/tables in file.
    
    Can act as context manger allowing use of 'with' statement.
    
    >>> with BinaryReader(binary_file) as reader:
            do something     
    ''' 
    def __init__(self, binary_file):
        self._file_handle = binary_file
    
    def __iter__(self):
        return self

    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        self.close()
    
    def next(self):
        '''
        Should return a table object.
        '''
        pass

class BinaryWriter(object):
    '''
    Class for writing BinaryFile object.
    
    Can act as context manger allowing use of 'with' statement.
    
    >>> with BinaryWriter(binary_file) as writer:
            do something
    
    ''' 
    def __init__(self, binary_file):
        self._file_handle = binary_file
    
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        self.close()
    
    def add_row(self, row):
        '''
        The format of row should be defined by inheriting class.
        
        All information about where the data is written should be figured out from the row object.
        '''
        pass
    
    def close(self):
        self._file_handle.close()

class BinaryTable(object):
    '''
    Key data object of a binary file. There does not need to be a one to one correspondence between tables in underlying
    HDF5 file and these objects.
    
    For example we could have the group chr1 and table chr1/index chr1/counts.
    '''
    def __init__(self, name):
        pass
    
    def __iter__(self):
        return self
    
    def next(self):
        '''
        Should be able to return a 'row' of the table.
        ''' 
        pass
    
class BinaryFile(object):
    '''
    Class representing a HDF5 file.
    
    Any access to the underlying HDF5 file hierarchy should be placed here.
    '''
    def __init__(self, file_name, file_mode, compression_level=1, compression_lib='zlib'):
        '''
        For compatibility it is recommended the compression values are left at defaults.
        
        Arguments:
        file_name -- Path to file
        file_mode -- How file should be opened i.e. r, w, a, r+
        compression_level -- Level of compression to use from 1 to 9
        compression_lib -- Compression library to use see PyTables docs for option.
        '''
        pass
