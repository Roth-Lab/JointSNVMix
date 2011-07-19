'''
Created on 2010-08-06

@author: Andrew Roth
'''
import random
import time

import numpy as np

from tables import openFile, Filters, UInt32Col, StringCol
from tables.description import IsDescription
    
class JointCountsReader(object):
    '''
    Class for reading jcnt files.
    '''
    def __init__(self, jcnt_file_name):
        self._file_handle = _JointCountsFile(jcnt_file_name, 'r')
    
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        self.close()
    
    def get_table_list(self):
        return self._file_handle.entries
    
    def get_counts(self, chrom=None):
        if chrom is None:
            counts = []
            
            for chrom in sorted(self.get_table_list()):
                counts.append(self._get_chrom_counts(chrom))
                
            counts = np.vstack(counts)
        else:
            counts = self._get_chrom_counts(chrom)

        
        return counts
    
    def get_random_counts_subsample(self, chrom, sample_size):
        table = self._file_handle.get_table(chrom)
        n = table.nrows
        
        table_sample_indices = random.sample(xrange(n), sample_size)
        
        sub_sample = table[table_sample_indices]
        
        return np.column_stack((
                              sub_sample['normal_counts_a'],
                              sub_sample['normal_counts_b'],
                              sub_sample['tumour_counts_a'],
                              sub_sample['tumour_counts_b']
                              ))
        
        
    
    def get_number_of_table_rows(self, chrom):
        table = self._file_handle.get_table(chrom)
        n = table.nrows
        
        return n
    
    def get_data_set_size(self):
        data_set_size = 0
        
        for chrom in self.get_table_list():
            data_set_size += self.get_number_of_table_rows(chrom)
            
        return data_set_size
    
    def get_table(self, chrom):
        return self._file_handle.get_table(chrom)
    
    def close(self):
        '''
        Should be called when done reading the file.
        '''
        self._file_handle.close()
    
    def _get_chrom_counts(self, chrom):
        table = self._file_handle.get_table(chrom)
            
        counts = np.column_stack([
                                    table.col('normal_counts_a'),
                                    table.col('normal_counts_b'),
                                    table.col('tumour_counts_a'),
                                    table.col('tumour_counts_b')
                                    ])
        
        return counts
    
class JointCountsWriter(object):
    '''
    Class for writing jcnt files.
    '''
    def __init__(self, jcnt_file_name):
        self._file_handle = _JointCountsFile(jcnt_file_name, 'w')
        
        self._row_buffer = {}
        self._buffer_size = 0
        
        # Max number of rows to buffer before writing to disk.
        self._max_buffer_size = 100000
    
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        self.close()
        
    def add_row(self, chrom, row):
        '''
        Add row to jcnt file under the chrom table. Row should be compatible with storage in _JointCountsIndexTable.
        
        For example
        
        [12345, 'A', 'C', 100, 0, 100, 100]
        '''        
        if chrom not in self._row_buffer:
            self._row_buffer[chrom] = []
        
        self._row_buffer[chrom].append(row)
        self._buffer_size += 1
        
        if self._buffer_size >= self._max_buffer_size:
            self._write_buffer()
            
    def close(self):
        '''
        Must be called when done if the class not used as context manager in with block.
        ''' 
        self._write_buffer()
        self._file_handle.close()
        
    def _write_buffer(self):
        for chrom in self._row_buffer:
            rows = self._row_buffer[chrom]
            self._file_handle.add_rows_to_table(chrom, rows)
            
        self._row_buffer = {}
        self._buffer_size = 0
    
class _JointCountsFile:
    '''
    Class representing a joint counts (jcnt) formated file. Should only be accessed using JointCountsReader and
    JointCountsWriter classes.
    
    Any access to the underlying HDF5 file hierachy should be placed here.
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
        compression_filters = Filters(complevel=compression_level, complib=compression_lib)

        if file_mode == "w":
            self._file_handle = openFile(file_name, file_mode, filters=compression_filters)
            
            self._file_handle.setNodeAttr('/', 'creation_date', time.ctime())
        else:
            self._file_handle = openFile(file_name, file_mode)

        self._init_entries()

        self._init_chrom_tables()
                   
    def add_rows_to_table(self, chrom, rows):
        '''
        Add rows to table chrom.
        '''
        table = self._get_chrom_table(chrom)
    
        table.append(rows)
        table.flush()

    def get_table(self, chrom):
        '''
        Get table chrom.
        '''
        return self._get_chrom_table(chrom)
        
    def close(self):
        self._file_handle.close()        

    def _get_chrom_table(self, chr_name):
        if chr_name in self._chrom_tables:
            chrom_table = self._chrom_tables[chr_name]
        else:
            chrom_table = self._file_handle.createTable(self._file_handle.root, chr_name, _JointCountsIndexTable)

            self._chrom_tables[chr_name] = chrom_table

        return chrom_table

    def _init_entries(self):
        entries = set()

        for table in self._file_handle.iterNodes(where=self._file_handle.root):
            entries.add(table._v_name)

        self.entries = entries

    def _init_chrom_tables(self):
        chrom_tables = {}

        for chr_name in self.entries:
            chrom_tables[chr_name] = self._file_handle.getNode(self._file_handle.root, chr_name)

        self._chrom_tables = chrom_tables

class _JointCountsIndexTable(IsDescription):
    position = UInt32Col(pos=0)

    ref_base = StringCol(itemsize=1, pos=1)

    non_ref_base = StringCol(itemsize=1, pos=2)

    normal_counts_a = UInt32Col(pos=3)
    
    normal_counts_b = UInt32Col(pos=4)
    
    tumour_counts_a = UInt32Col(pos=5)
    
    tumour_counts_b = UInt32Col(pos=6)
