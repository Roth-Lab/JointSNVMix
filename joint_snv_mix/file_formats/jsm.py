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
        
    def get_rows(self, chr_name):
        return self._file_handle.get_rows(chr_name)
    
    def get_parameters(self):
        return self._file_handle.get_parameters()        
               
class JointSnvMixWriter:
    def __init__(self, file_name,):
        self._file_handle = _JointSnvMixFile(file_name, 'w')
        
    def write_priors(self, priors):
        self._file_handle.write_priors(priors)
        
    def write_parameters(self, parameters):
        self._file_handle.write_parameters(parameters)
        
    def add_row(self, chr_name, jcnt_row):
        data = []
        
        for jcnt_row, resp in zip(jcnt_rows.tolist(), responsibilities):
            row = []
            row.extend(jcnt_row)
            row.extend(resp)
            data.append(row)
        
        self._file_handle.write_chr_table(chr_name, data)

    def close(self):
        self._file_handle.close()
        
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
                
    def get_priors(self):
        priors = {}
        
        self._read_tree(priors, self._priors_group)
        
        return priors
    
    def get_parameters(self):
        parameters = {}
        
        self._read_tree(parameters, self._parameters_group)
        
        return parameters
    
    def _read_tree(self, params, group):
        for entry in self._file_handle.iterNodes(where=group):
            name = entry._v_name 
            
            if isinstance(entry, Leaf):
                params[name] = entry[:]
            else:
                params[name] = {}
                self._read_tree(params[name], entry)

    def write_chr_table(self, chr_name, data):
        if chr_name not in self._chrom_tables:
            chr_table = self._file_handle.createTable('/data', chr_name, JointSnvMixTable)
            
            self._chrom_tables[chr_name] = chr_table
        else:
            chr_table = self._chrom_tables[chr_name]        
        
        chr_table.append(data)
        
    def get_responsibilities(self, chr_name):
        table = self._chrom_tables[chr_name]
        
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
    
    def get_rows(self, chr_name, row_indices=None):
        table = self._chrom_tables[chr_name]
        
        if row_indices is None:
            return table
        else:
            return table[row_indices]
    
    def get_position(self, chr_name, coord):
        table = self._chrom_tables[chr_name]
        
        search_string = "position == {0}".format(coord)
        row = table.readWhere(search_string)
        
        if len(row) == 0:
            row = []
        else:
            row = row[0].tolist()
        
        return row
        
        
    def close(self):
        self._file_handle.close()

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
    
class JointSnvMixTable(IsDescription):
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
