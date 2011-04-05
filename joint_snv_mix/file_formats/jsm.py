'''
Created on 2010-08-25

@author: Andrew Roth
'''
import time

import numpy as np

from tables import openFile, Filters, Float64Atom, StringCol, IsDescription, UInt32Col, Float64Col, Leaf
   
class JointSnvMixReader:
    def __init__(self, file_name):
        self._file_handle = _JointSnvMixFile(file_name, 'r')

    def get_table_list(self):
        return self._file_handle.entries
    
    def close(self):
        self._file_handle.close()
        
    def get_table(self, chr_name):
        return self._file_handle.get_table(chr_name)
    
    def get_parameters(self):
        return self._file_handle.get_parameters()        
               
class JointSnvMixWriter:
    '''
    Class for writing jcnt files.
    '''
    def __init__(self, jcnt_file_name):
        self._file_handle = _JointSnvMixFile(jcnt_file_name, 'w')
        
        self._row_buffer = {}
        self._buffer_size = 0
        
        # Max number of rows to buffer before writing to disk.
        self._max_buffer_size = 100000
    
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        self.close()
        
    def write_priors(self, priors):
        self._file_handle.write_priors(priors)
        
    def write_parameters(self, parameters):
        self._file_handle.write_parameters(parameters)
        
    def add_row(self, chrom, row):
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
        
class _JointSnvMixFile:
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
            self._data_group = self._file_handle.createGroup("/", "data")
            self._parameters_group = self._file_handle.createGroup("/", "parameters")
            self._priors_group = self._file_handle.createGroup("/", "priors")
            
            self._file_handle.setNodeAttr('/', 'creation_date', time.ctime())
        else:
            self._file_handle = openFile(file_name, file_mode)
            self._data_group = self._file_handle.root.data
            self._parameters_group = self._file_handle.root.parameters
            self._priors_group = self._file_handle.root.priors

        self._init_entries()

        self._init_chr_tables()
        
    def write_priors(self, priors):
        priors_group = self._priors_group
        
        self._write_tree(priors, priors_group)
        
    def write_parameters(self, parameters):    
        params_group = self._parameters_group
        
        self._write_tree(parameters, params_group)
                
    def get_priors(self):
        priors = {}
        
        self._read_tree(priors, self._priors_group)
        
        return priors
    
    def get_parameters(self):
        parameters = {}
        
        self._read_tree(parameters, self._parameters_group)
        
        return parameters

    def add_rows_to_table(self, chrom, rows):
        '''
        Add rows to table chrom.
        '''
        table = self._get_chrom_table(chrom)
    
        table.append(rows)
        table.flush()
        
    def get_responsibilities(self, chrom):
        table = self._chrom_tables[chrom]
        
        responsibilities = np.column_stack((
                                            table.col('p_aa_aa'),
                                            table.col('p_aa_ab'),
                                            table.col('p_aa_bb'),
                                            table.col('p_ab_aa'),
                                            table.col('p_ab_ab'),
                                            table.col('p_ab_bb'),
                                            table.col('p_bb_aa'),
                                            table.col('p_bb_ab'),
                                            table.col('p_bb_bb')
                                            ))
        
        return responsibilities
    
    def get_table(self, chrom, row_indices=None):
        table = self._chrom_tables[chrom]
        
        if row_indices is None:
            return table
        else:
            return table[row_indices]
    
    def get_position(self, chrom, coord):
        table = self._chrom_tables[chrom]
        
        search_string = "position == {0}".format(coord)
        row = table.readWhere(search_string)
        
        if len(row) == 0:
            row = []
        else:
            row = row[0].tolist()
        
        return row
        
        
    def close(self):
        self._file_handle.close()

    def _write_tree(self, params, group):
        for name, value in params.items():
            if isinstance(value, dict):
                new_group = self._file_handle.createGroup(group, name)
                self._write_tree(value, new_group)
            else:
                atom = Float64Atom(())
        
                shape = np.array(value).shape
        
                parameter_array = self._file_handle.createCArray(group, name, atom, shape)
                
                parameter_array[:] = value[:]
    
    def _read_tree(self, params, group):
        for entry in self._file_handle.iterNodes(where=group):
            name = entry._v_name 
            
            if isinstance(entry, Leaf):
                params[name] = entry[:]
            else:
                params[name] = {}
                self._read_tree(params[name], entry)
                
    
    def _get_chrom_table(self, chrom):
        if chrom in self._chrom_tables:
            chrom_table = self._chrom_tables[chrom]
        else:
            chrom_table = self._file_handle.createTable(self._data_group, chrom, _JointSnvMixTable)

            self._chrom_tables[chrom] = chrom_table

        return chrom_table

    def _init_entries(self):
        '''
        Build the initial list of chromosomes in table.
        '''
        entries = set()

        for table in self._file_handle.iterNodes(where=self._data_group):
            entries.add(table._v_name)

        self.entries = entries

    def _init_chr_tables(self):
        chr_tables = {}

        for chr_name in self.entries:
            chr_tables[chr_name] = self._file_handle.getNode(self._data_group, chr_name)

        self._chrom_tables = chr_tables
    
class _JointSnvMixTable(IsDescription):
    position = UInt32Col(pos=0)

    ref_base = StringCol(itemsize=1, pos=1)

    normal_base = StringCol(itemsize=1, pos=2)

    tumour_base = StringCol(itemsize=1, pos=3)

    normal_counts_a = UInt32Col(pos=4)
    
    normal_counts_b = UInt32Col(pos=5)
    
    tumour_counts_a = UInt32Col(pos=6)
    
    tumour_counts_b = UInt32Col(pos=7)
    
    p_aa_aa = Float64Col(pos=8)
    
    p_aa_ab = Float64Col(pos=9)
    
    p_aa_bb = Float64Col(pos=10)
    
    p_ab_aa = Float64Col(pos=11)
    
    p_ab_ab = Float64Col(pos=12)
    
    p_ab_bb = Float64Col(pos=13)
    
    p_bb_aa = Float64Col(pos=14)
    
    p_bb_ab = Float64Col(pos=15)
    
    p_bb_bb = Float64Col(pos=16) 
