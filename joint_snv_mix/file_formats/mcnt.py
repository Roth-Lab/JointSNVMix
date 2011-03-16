'''
Created on 2010-08-06

@author: Andrew Roth
'''
import time

import numpy as np

from tables import openFile, Filters, UInt32Col, StringCol
from tables.description import IsDescription

class MultinomialCountsFile:
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

        self._chrom_tables = self._init_chr_tables()
    
    def add_rows( self, chr_name, rows ):
        table = self._get_chr_table( chr_name )
        
        table.append( rows )
        table.flush()
        
    def get_rows( self, chr_name, indices=None ):
        table = self._get_chr_table( chr_name )
        
        if indices is None:
            rows = table.read()
        else:
            rows = table[indices]
        
        return rows
    
    def get_table_size( self, chr_name ):
        table = self._get_chr_table( chr_name )
        
        return table.nrows
    
    def close( self ):
        self._file_handle.close()

    def _get_chr_table( self, chr_name ):
        '''
        Get a counts table for a chromosome.
        
        Fetch the table if it exists otherwise create it.
        
        Arguments:
        chr_name -- Name of table to fetch.
        
        Return:
        chr_table -- A counts table object. See BinomialCountsTable for columns.
        '''
        if chr_name in self._chrom_tables:
            chr_table = self._chrom_tables[chr_name]
        else:
            chr_table = self._file_handle.createTable( self._file_handle.root, chr_name, MultinomialCountsIndexTable )

            self._chrom_tables[chr_name] = chr_table

        return chr_table

    def _init_entries( self ):
        '''
        Build the initial list of chromosomes in table.
        '''
        entries = set()

        for table in self._file_handle.iterNodes( where=self._file_handle.root ):
            entries.add( table._v_name )

        return entries

    def _init_chr_tables( self ):
        chr_tables = {}

        for chr_name in self.entries:
            chr_tables[chr_name] = self._file_handle.getNode( self._file_handle.root, chr_name )

        return chr_tables

class MultinomialCountsReader:
    '''
    Helper class to simpilfy reading jcnt files.
    '''
    def __init__( self, file_name ):
        '''
        Arguments:
        file_name -- Path to joint counts file to be read.
        '''
        self._file_handle = MultinomialCountsFile( file_name, 'r' )
        
    def close( self ):
        '''
        Should be called when reading the file.
        '''
        self._file_handle.close()
    
    def get_chr_list( self ):
        return self._file_handle.entries
    
    def get_counts( self, chr_name=None, indices=None ):
        if chr_name is None:
            counts = []
            
            for chr_name in sorted( self.get_chr_list() ):
                rows = self._file_handle.get_rows( chr_name )
                
                counts_cols = np.column_stack( [
                                         rows['normal_counts_A'],
                                         rows['normal_counts_C'],
                                         rows['normal_counts_G'],
                                         rows['normal_counts_T'],
                                         rows['tumour_counts_A'],
                                         rows['tumour_counts_C'],
                                         rows['tumour_counts_G'],
                                         rows['tumour_counts_T']
                                         ] )
                
                counts.append( counts_cols )
                
            counts = np.vstack( counts )
        else:
            if indices is None:
                rows = self._file_handle.get_rows( chr_name )
            else:
                rows = self._file_handle.get_rows( chr_name, indices )
                
            counts = np.column_stack( [
                             rows['normal_counts_A'],
                             rows['normal_counts_C'],
                             rows['normal_counts_G'],
                             rows['normal_counts_T'],
                             rows['tumour_counts_A'],
                             rows['tumour_counts_C'],
                             rows['tumour_counts_G'],
                             rows['tumour_counts_T']
                        ] )
        
        return counts
    
    def get_chr_size( self, chr_name ):
        return self._file_handle.get_table_size( chr_name )
    
    def get_data_set_size( self ):
        data_set_size = 0
        
        for chr_name in self.get_chr_list():
            data_set_size += self._file_handle.get_table_size( chr_name )
            
        return data_set_size
    
    def get_rows( self, chr_name, indices=None ):
        return self._file_handle.get_rows( chr_name, indices )

class MultinomialCountsIndexTable( IsDescription ):
    position = UInt32Col( pos=0 )

    ref_base = StringCol( itemsize=1, pos=1 )

    normal_counts_A = UInt32Col( pos=2 )
    
    normal_counts_C = UInt32Col( pos=3 )
    
    normal_counts_G = UInt32Col( pos=4 )
    
    normal_counts_T = UInt32Col( pos=5 )
    
    tumour_counts_A = UInt32Col( pos=6 )
    
    tumour_counts_C = UInt32Col( pos=7 )
    
    tumour_counts_G = UInt32Col( pos=8 )
    
    tumour_counts_T = UInt32Col( pos=9 )
