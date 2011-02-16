'''
Created on 2010-08-06

@author: Andrew Roth
'''
import time

import numpy as np

from tables import openFile, Filters, UInt32Col, StringCol
from tables.description import IsDescription

class ConanCountsFile:
    '''
    Class representing a joint counts formated file.
    
    Any acess to the underlying HDF5 file hierachy should be placed here.
    '''
    def __init__( self, file_name, file_mode, compression_level=1, compression_lib='zlib' ):
        '''
        For compatibility it is recommended the compression values are left at defaults.
        
        Arguments:
        file_name -- Path to file
        file_mode -- How file should be opened i.e. r, w, a, r+
        compression_level -- Level of compression to use from 1 to 9
        compression_lib -- Compression library to use see PyTables docs for option.
        '''
        compression_filters = Filters( complevel=compression_level, complib=compression_lib )

        if file_mode == "w":
            self._file_handle = openFile( file_name, file_mode, filters=compression_filters )
            
            self._file_handle.setNodeAttr( '/', 'creation_date', time.ctime() )
        else:
            self._file_handle = openFile( file_name, file_mode )

        self.entries = self._init_entries()

        self._cn_groups = self._init_cn_groups()
    
    def add_rows( self, cn_status, chr_name, rows ):
        table = self._get_chr_table( cn_status, chr_name )
        
        table.append( rows )
        table.flush()
        
    def get_rows( self, cn_status, chr_name ):
        table = self._get_chr_table( cn_status, chr_name )
        
        rows = table.read()
        
        return rows
    
    def get_table_size( self, cn_status, chr_name ):
        table = self._get_chr_table( cn_status, chr_name )
        
        return table.nrows
    
    def get_chr_list( self, cn_status ):
        return self.entries[cn_status]
    
    def close( self ):
        self._file_handle.close()

    def _get_chr_table( self, cn_status, chr_name ):
        '''
        Get a counts table for a chromosome.
        
        Fetch the table if it exists otherwise create it.
        
        Arguments:
        chr_name -- Name of table to fetch.
        
        Return:
        chr_table -- A counts table object. See BinomialCountsTable for columns.
        '''
        if cn_status in self._cn_groups:
            cn_group = self._cn_groups[cn_status]
        else:
            cn_group = self._file_handle.createGroup( self._file_handle.root, cn_status )

            self._cn_groups[cn_status] = cn_group
        
        if chr_name in cn_group:
            chr_table = self._file_handle.getNode( cn_group, chr_name )
        else:
            chr_table = self._file_handle.createTable( cn_group, chr_name, JointCountsIndexTable )

        return chr_table

    def _init_entries( self ):
        '''
        Build the initial list of chromosomes in table.
        '''
        entries = {}

        for cn_group in self._file_handle.iterNodes( where=self._file_handle.root ):
            group_name = cn_group._v_name
            
            cn_tables = self._file_handle.listNodes( where=cn_group )
            entries[group_name] = set( [x._v_name for x in cn_tables] )

        return entries

    def _init_cn_groups( self ):
        cn_groups = {}

        for cn_state in self.entries:
            cn_groups[cn_state] = self._file_handle.getNode( self._file_handle.root, cn_state )

        return cn_groups

class ConanCountsReader:
    '''
    Helper class to simpilfy reading jcnt files.
    '''
    def __init__( self, file_name ):
        '''
        Arguments:
        file_name -- Path to joint counts file to be read.
        '''
        self._file_handle = ConanCountsFile( file_name, 'r' )
        
    def close( self ):
        '''
        Should be called when reading the file.
        '''
        self._file_handle.close()
    
    def get_cn_states( self ):
        return self._file_handle.entries.keys()
    
    def get_chr_list( self, cn_state ):
        return self._file_handle.get_chr_list( cn_state )
    
    def get_counts( self, cn_state=None, chr_name=None ):
        if cn_state is None:
            counts = self._load_counts()
        elif chr_name is None:
            counts = self._load_cn_counts( cn_state )            
        else:
            counts = self._load_chr_counts( cn_state, chr_name )
        
        return counts
    
    def _load_chr_counts( self, cn_state, chr_name ):
        rows = self._file_handle.get_rows( cn_state, chr_name )
        
        counts_cols = np.column_stack( [
                                 rows['normal_counts_a'], rows['normal_counts_b'],
                                 rows['tumour_counts_a'], rows['tumour_counts_b']
                                 ] )
        
        return counts_cols
    
    def _load_cn_counts( self, cn_state ):
        counts = []
        
        chr_list = self._file_handle.get_chr_list( cn_state )
        
        for chr_name in sorted( chr_list ):                    
            counts_cols = self._load_chr_counts( cn_state, chr_name )
            
            counts.append( counts_cols )
        
        counts = np.vstack( counts )
        
        return counts
    
    def _load_counts( self ):
        counts = []
        
        cn_states = self.get_cn_states()
        
        for cn_state in sorted( cn_states ):
            counts_cols = self._load_cn_counts( cn_state )

            counts.append( counts_cols )
        
        counts = np.vstack( counts )
        
        return counts
    
    def get_chr_size( self, cn_status, chr_name ):
        return self._file_handle.get_table_size( cn_status, chr_name )
    
    def get_data_set_size( self, cn_state=None ):
        data_set_size = 0
        
        if cn_state is None:        
            for cn_state in self.get_cn_states():
                for chr_name in self._file_handle.get_chr_list( cn_state ):
                    data_set_size += self._file_handle.get_table_size( cn_state, chr_name )
        else:
            for chr_name in self._file_handle.get_chr_list( cn_state ):
                    data_set_size += self._file_handle.get_table_size( cn_state, chr_name )
            
        return data_set_size
        
    def get_rows( self, cn_state, chr_name ):
        return self._file_handle.get_rows( cn_state, chr_name )

class JointCountsIndexTable( IsDescription ):
    position = UInt32Col( pos=0 )

    ref_base = StringCol( itemsize=1, pos=1 )

    normal_base = StringCol( itemsize=1, pos=2 )

    tumour_base = StringCol( itemsize=1, pos=3 )

    normal_counts_a = UInt32Col( pos=4 )
    
    normal_counts_b = UInt32Col( pos=5 )
    
    tumour_counts_a = UInt32Col( pos=6 )
    
    tumour_counts_b = UInt32Col( pos=7 )
