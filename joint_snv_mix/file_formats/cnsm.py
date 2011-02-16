'''
Created on 2010-08-25

@author: Andrew Roth
'''
import time

from joint_snv_mix.file_formats.jsm import JointSnvMixTable
from tables import Float64Atom, StringCol, UInt32Col, Float64Col, Filters, openFile
import numpy as np
from joint_snv_mix.file_formats.cncnt import JointCountsIndexTable
from joint_snv_mix import constants

   
class ConanSnvMixFile:
    def __init__( self, file_name, file_mode, compression_level=1, compression_lib='zlib' ):
        '''Constructor
                    
        For compatibility it is reccomeded the compression values are left at defaults.
        
        Arguments:
        file_name -- Path to file
        file_mode -- How file should be opened i.e. r, w, a, r+
        compression_level -- Level of compression to use from 1 to 9
        compression_lib -- Compression library to use see PyTables docs for option.
        '''
        compression_filters = Filters( complevel=compression_level, complib=compression_lib )

        if file_mode == "w":
            self._file_handle = openFile( file_name, file_mode, filters=compression_filters )
            self._data_group = self._file_handle.createGroup( "/", "data" )
            self._parameters_group = self._file_handle.createGroup( "/", "parameters" )
            self._priors_group = self._file_handle.createGroup( "/", "priors" )
            
            self._file_handle.setNodeAttr( '/', 'creation_date', time.ctime() )
        else:
            self._file_handle = openFile( file_name, file_mode )
            self._data_group = self._file_handle.root.data
            self._parameters_group = self._file_handle.root.parameters
            self._priors_group = self._file_handle.root.priors

        self._init_entries()

        self._init_cn_groups()
    
    def get_chr_list( self, cn_status ):
        return self.entries[cn_status]
    
    def write_priors( self, priors ):
        priors_group = self._priors_group
        
        self._write_tree( priors, priors_group )
        
    def write_parameters( self, parameters ):    
        params_group = self._parameters_group
        
        self._write_tree( parameters, params_group )
    
    def _write_tree( self, params, group ):
        for name, value in params.items():
            if isinstance( value, dict ):
                new_group = self._file_handle.createGroup( group, name )
                self._write_tree( value, new_group )
            else:
                atom = Float64Atom( () )
        
                shape = np.array( value ).shape
        
                parameter_array = self._file_handle.createCArray( group, name, atom, shape )
                
                parameter_array[:] = value[:]
                
    def get_priors( self ):
        priors = {}
        
        self._read_tree( priors, self._priors_group )
        
        return priors
    
    def get_parameters( self ):
        parameters = {}
        
        self._read_tree( parameters, self._parameters_group )
        
        return parameters
        
    def write_chr_table( self, cn_state, chr_name, index_rows, soft_labels ):
        if cn_state in self._cn_groups:
            cn_group = self._cn_groups[cn_state]
        else:
            cn_group = self._file_handle.createGroup( '/data', cn_state )

            self._cn_groups[cn_state] = cn_group
        
        if chr_name in cn_group:
            chr_group = self._file_handle.getNode( cn_group, chr_name )
        else:
            chr_group = self._file_handle.createGroup( cn_group, chr_name )
            
            self._file_handle.createTable( chr_group, 'index' , JointCountsIndexTable )
            
            atom = Float64Atom( () )
            nclass = 3 * constants.cn_state_map[cn_state]
            shape = ( 0, nclass )
            
            self._file_handle.createEArray( chr_group, 'soft_labels', atom, shape )

        
        chr_group.index.append( index_rows )
        chr_group.soft_labels.append( soft_labels )
        
    def get_responsibilities( self, cn_state, chr_name ):
        cn_group = self._cn_groups[cn_state]
        table = self._file_handle.getnode( cn_group, chr_name )
        
        responsibilities = np.column_stack( ( 
                                            table.col( 'p_aa_aa' ),
                                            table.col( 'p_aa_ab' ),
                                            table.col( 'p_aa_bb' ),
                                            table.col( 'p_ab_aa' ),
                                            table.col( 'p_ab_ab' ),
                                            table.col( 'p_ab_bb' ),
                                            table.col( 'p_bb_aa' ),
                                            table.col( 'p_bb_ab' ),
                                            table.col( 'p_bb_bb' )
                                            ) )
        
        return responsibilities
    
    def get_rows( self, cn_state, chr_name, row_indices=None ):
        cn_group = self._cn_groups[cn_state]
        chr_group = self._file_handle.getNode( cn_group, chr_name )
        
        if row_indices is None:
            return chr_group.index, chr_group.soft_labels
        else:
            return chr_group.index[row_indices], chr_group.soft_labels[row_indices]
    
    def get_position( self, chr_name, coord ):
        '''
        Search file for a given position.
        '''
        for cn_state in self._cn_groups:            
            cn_group = self._cn_groups[cn_state]
            table = self._file_handle.getnode( cn_group, chr_name )
            
            search_string = "position == {0}".format( coord )
            row = table.readWhere( search_string )
            
            if len( row ) == 0:
                row = []
            else:
                row = row[0].tolist()
                row.append( cn_state )
                return row
            
            return row
                
    def close( self ):
        self._file_handle.close()

    def _init_entries( self ):
        '''
        Build the initial list of chromosomes in table.
        '''
        entries = {}

        for cn_group in self._file_handle.iterNodes( where='/data' ):
            group_name = cn_group._v_name
            
            cn_tables = self._file_handle.listNodes( where=cn_group )
            entries[group_name] = set( [x._v_name for x in cn_tables] )

        self.entries = entries

    def _init_cn_groups( self ):
        cn_groups = {}

        for cn_state in self.entries:
            cn_groups[cn_state] = self._file_handle.getNode( '/data', cn_state )

        self._cn_groups = cn_groups


class ConanSnvMixReader:
    def __init__( self, file_name ):
        self._file_handle = ConanSnvMixFile( file_name, 'r' )

    def get_cn_states( self ):
        return self._file_handle.entries.keys()
    
    def get_chr_list( self, cn_state ):
        return self._file_handle.get_chr_list( cn_state )
        
    def get_position( self, chr_name, coord ):
        return self._file_handle.get_position( chr_name, coord )
    
    def get_parameters( self ):
        return self._file_handle.get_parameters()     

    def get_rows( self, cn_state, chr_name ):
        return self._file_handle.get_rows( cn_state, chr_name )

    def close( self ):
        self._file_handle.close()
               
class ConanSnvMixWriter:
    def __init__( self, file_name, ):
        self._file_handle = ConanSnvMixFile( file_name, 'w' )
        
    def write_priors( self, priors ):
        self._file_handle.write_priors( priors )
        
    def write_parameters( self, parameters ):
        self._file_handle.write_parameters( parameters )
        
    def write_data( self, cn_state, chr_name, jcnt_rows, responsibilities ):        
        self._file_handle.write_chr_table( cn_state, chr_name, jcnt_rows, responsibilities )

    def close( self ):
        self._file_handle.close()
