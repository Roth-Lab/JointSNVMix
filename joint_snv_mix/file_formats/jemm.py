'''
Created on 2010-08-25

@author: Andrew Roth
'''
import time

import numpy as np

from tables import openFile, Filters, Float64Atom, StringCol, IsDescription, UInt32Col, Float64Col
from joint_snv_mix.constants import joint_extended_multinomial_genotypes
import joint_snv_mix.constants as constants
   
class JointExtendedMultiMixFile:
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

        self._init_chr_tables()
        
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
    
    def _read_tree( self, params, group ):
        for entry in self._file_handle.iterNodes( where=group ):
            name = entry._v_name 
            
            if isinstance( entry, Leaf ):
                params[name] = entry[:]
            else:
                params[name] = {}
                self._read_tree( params[name], entry )
            
    def write_chr_table( self, chr_name, data ):
        if chr_name not in self._chrom_tables:
            chr_table = self._file_handle.createTable( '/data', chr_name, JointSnvMixTable )
            
            self._chrom_tables[chr_name] = chr_table
        else:
            chr_table = self._chrom_tables[chr_name]
        
        chr_table.append( data )
        
    def get_responsibilities( self, chr_name ):
        table = self._chrom_tables[chr_name]
        
        probs = []
        
        for jmg in joint_extended_multinomial_genotypes:
            prob = "_".join( jmg )
            prob = "p" + "_" + prob
            
            probs.append( table.col( prob ) )
        
        responsibilities = np.column_stack( probs )
        
        return responsibilities
    
    def get_row_above_prob( self, chr_name, class_labels, prob_threshold ):
        table = self._chrom_tables[chr_name]
        
        probs = []
        
        for i in class_labels:
            prob = "_".join( joint_extended_multinomial_genotypes[i] )
            prob = "p" + "_" + prob
            
            probs.append( prob )
            
        
        query_string = " + ".join( probs )
        
        query_string = "{0} >= {1}".format( query_string, prob_threshold )
        
        rows = table.readWhere( query_string )
        
        return rows
    
    def get_rows( self, chr_name, row_indices=None ):
        table = self._chrom_tables[chr_name]
        
        if row_indices is None:
            return table[:]
        else:
            return table[row_indices]
    
    def get_position( self, chr_name, coord ):
        table = self._chrom_tables[chr_name]
        
        search_string = "position == {0}".format( coord )
        row = table.readWhere( search_string )
        
        if len( row ) == 0:
            row = []
        else:
            row = row[0].tolist()
        
        return row
        
        
    def close( self ):
        self._file_handle.close()

    def _init_entries( self ):
        '''
        Build the initial list of chromosomes in table.
        '''
        entries = set()

        for table in self._file_handle.iterNodes( where=self._data_group ):
            entries.add( table._v_name )

        self.entries = entries

    def _init_chr_tables( self ):
        chr_tables = {}

        for chr_name in self.entries:
            chr_tables[chr_name] = self._file_handle.getNode( self._data_group, chr_name )

        self._chrom_tables = chr_tables

class JointExtendedMultiMixReader:
    def __init__( self, file_name ):
        self._file_handle = JointExtendedMultiMixFile( file_name, 'r' )

    def get_chr_list( self ):
        return self._file_handle.entries
    
    def get_genotype_rows_by_argmax( self, chr_name, genotype_class ):
        if genotype_class == 'Somatic':
            class_labels = constants.somatic_extended_multinomial_genotypes_indices
        elif genotype_class == 'Germline':
            class_labels = constants.matched_extended_multinomial_genotypes_indices
        elif genotype_class == 'LOH':
            class_labels = constants.loh_extended_multinomial_genotypes_indices
        else:
            raise Exception( 'Class {0} not accepted.'.format( genotype_class ) )
        
        rows = self._get_rows_by_argmax( chr_name, class_labels )
        
        return rows
    
    def get_genotype_rows_by_prob( self, chr_name, genotype_class, prob_threshold ):
        if genotype_class == 'Somatic':
            class_labels = constants.somatic_extended_multinomial_genotypes_indices
        elif genotype_class == 'Germline':
            class_labels = constants.matched_extended_multinomial_genotypes_indices
        elif genotype_class == 'LOH':
            class_labels = constants.loh_extended_multinomial_genotypes_indices
        else:
            raise Exception( 'Class {0} not accepted.'.format( genotype_class ) )
        
        rows = self._file_handle.get_row_above_prob( chr_name, class_labels, prob_threshold )

        return rows
    
    def get_position( self, chr_name, coord ):
        return self._file_handle.get_position( chr_name, coord )
    
    def close( self ):
        self._file_handle.close()
        
    def get_rows( self, chr_name ):
        return self._file_handle.get_rows( chr_name )
    
    def _get_rows_by_argmax( self, chr_name, class_labels ):
        responsibilities = self._file_handle.get_responsibilities( chr_name )
        
        labels = np.argmax( responsibilities, axis=1 )
        
        row_indices = []
        for class_label in class_labels:
            row_indices.extend( np.where( labels == class_label )[0] )
        
        row_indices = sorted( row_indices )
        
        rows = self._file_handle.get_rows( chr_name, row_indices )
        
        return rows
    
    def _get_rows_by_prob( self, chr_name, class_labels, prob_threshold ):
        responsibilities = self._file_handle.get_responsibilities( chr_name )
        
        shape = ( responsibilities.shape[0], )
        class_prob = np.zeros( shape )
        
        for class_label in class_labels:
            class_prob += responsibilities[:, class_label]
        
        row_indices = np.where( class_prob >= prob_threshold )
        
        if len( row_indices ) > 0:
            row_indices = row_indices[0]
            rows = self._file_handle.get_rows( chr_name, row_indices )
            
            return rows
        else:
            return []
               
class JointExtendedMultiMixWriter:
    def __init__( self, file_name, ):
        self._file_handle = JointExtendedMultiMixFile( file_name, 'w' )
        
    def write_priors( self, priors ):
        self._file_handle.write_priors( priors )
        
    def write_parameters( self, parameters ):
        self._file_handle.write_parameters( parameters )
        
    def write_data( self, chr_name, jcnt_rows, responsibilities ):
        data = []
        
        for jcnt_row, resp in zip( jcnt_rows.tolist(), responsibilities ):
            row = []
            row.extend( jcnt_row )
            row.extend( resp )
            data.append( row )
        
        self._file_handle.write_chr_table( chr_name, data )

    def close( self ):
        self._file_handle.close()
    
class JointExtendedMultiMixTable( IsDescription ):
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
    
    p_AA_AA = Float64Col( pos=10 )
    p_AA_AC = Float64Col( pos=11 )
    p_AA_AG = Float64Col( pos=12 )
    p_AA_AT = Float64Col( pos=13 )
    p_AA_CC = Float64Col( pos=14 )
    p_AA_CG = Float64Col( pos=15 )
    p_AA_CT = Float64Col( pos=16 )
    p_AA_GG = Float64Col( pos=17 )
    p_AA_GT = Float64Col( pos=18 )
    p_AA_TT = Float64Col( pos=19 )
    p_AA_ACG = Float64Col( pos=20 )
    p_AA_ACT = Float64Col( pos=21 )
    p_AA_AGT = Float64Col( pos=22 )
    p_AA_CGT = Float64Col( pos=23 )
    p_AA_ACGT = Float64Col( pos=24 )
    p_AC_AA = Float64Col( pos=25 )
    p_AC_AC = Float64Col( pos=26 )
    p_AC_AG = Float64Col( pos=27 )
    p_AC_AT = Float64Col( pos=28 )
    p_AC_CC = Float64Col( pos=29 )
    p_AC_CG = Float64Col( pos=30 )
    p_AC_CT = Float64Col( pos=31 )
    p_AC_GG = Float64Col( pos=32 )
    p_AC_GT = Float64Col( pos=33 )
    p_AC_TT = Float64Col( pos=34 )
    p_AC_ACG = Float64Col( pos=35 )
    p_AC_ACT = Float64Col( pos=36 )
    p_AC_AGT = Float64Col( pos=37 )
    p_AC_CGT = Float64Col( pos=38 )
    p_AC_ACGT = Float64Col( pos=39 )
    p_AG_AA = Float64Col( pos=40 )
    p_AG_AC = Float64Col( pos=41 )
    p_AG_AG = Float64Col( pos=42 )
    p_AG_AT = Float64Col( pos=43 )
    p_AG_CC = Float64Col( pos=44 )
    p_AG_CG = Float64Col( pos=45 )
    p_AG_CT = Float64Col( pos=46 )
    p_AG_GG = Float64Col( pos=47 )
    p_AG_GT = Float64Col( pos=48 )
    p_AG_TT = Float64Col( pos=49 )
    p_AG_ACG = Float64Col( pos=50 )
    p_AG_ACT = Float64Col( pos=51 )
    p_AG_AGT = Float64Col( pos=52 )
    p_AG_CGT = Float64Col( pos=53 )
    p_AG_ACGT = Float64Col( pos=54 )
    p_AT_AA = Float64Col( pos=55 )
    p_AT_AC = Float64Col( pos=56 )
    p_AT_AG = Float64Col( pos=57 )
    p_AT_AT = Float64Col( pos=58 )
    p_AT_CC = Float64Col( pos=59 )
    p_AT_CG = Float64Col( pos=60 )
    p_AT_CT = Float64Col( pos=61 )
    p_AT_GG = Float64Col( pos=62 )
    p_AT_GT = Float64Col( pos=63 )
    p_AT_TT = Float64Col( pos=64 )
    p_AT_ACG = Float64Col( pos=65 )
    p_AT_ACT = Float64Col( pos=66 )
    p_AT_AGT = Float64Col( pos=67 )
    p_AT_CGT = Float64Col( pos=68 )
    p_AT_ACGT = Float64Col( pos=69 )
    p_CC_AA = Float64Col( pos=70 )
    p_CC_AC = Float64Col( pos=71 )
    p_CC_AG = Float64Col( pos=72 )
    p_CC_AT = Float64Col( pos=73 )
    p_CC_CC = Float64Col( pos=74 )
    p_CC_CG = Float64Col( pos=75 )
    p_CC_CT = Float64Col( pos=76 )
    p_CC_GG = Float64Col( pos=77 )
    p_CC_GT = Float64Col( pos=78 )
    p_CC_TT = Float64Col( pos=79 )
    p_CC_ACG = Float64Col( pos=80 )
    p_CC_ACT = Float64Col( pos=81 )
    p_CC_AGT = Float64Col( pos=82 )
    p_CC_CGT = Float64Col( pos=83 )
    p_CC_ACGT = Float64Col( pos=84 )
    p_CG_AA = Float64Col( pos=85 )
    p_CG_AC = Float64Col( pos=86 )
    p_CG_AG = Float64Col( pos=87 )
    p_CG_AT = Float64Col( pos=88 )
    p_CG_CC = Float64Col( pos=89 )
    p_CG_CG = Float64Col( pos=90 )
    p_CG_CT = Float64Col( pos=91 )
    p_CG_GG = Float64Col( pos=92 )
    p_CG_GT = Float64Col( pos=93 )
    p_CG_TT = Float64Col( pos=94 )
    p_CG_ACG = Float64Col( pos=95 )
    p_CG_ACT = Float64Col( pos=96 )
    p_CG_AGT = Float64Col( pos=97 )
    p_CG_CGT = Float64Col( pos=98 )
    p_CG_ACGT = Float64Col( pos=99 )
    p_CT_AA = Float64Col( pos=100 )
    p_CT_AC = Float64Col( pos=101 )
    p_CT_AG = Float64Col( pos=102 )
    p_CT_AT = Float64Col( pos=103 )
    p_CT_CC = Float64Col( pos=104 )
    p_CT_CG = Float64Col( pos=105 )
    p_CT_CT = Float64Col( pos=106 )
    p_CT_GG = Float64Col( pos=107 )
    p_CT_GT = Float64Col( pos=108 )
    p_CT_TT = Float64Col( pos=109 )
    p_CT_ACG = Float64Col( pos=110 )
    p_CT_ACT = Float64Col( pos=111 )
    p_CT_AGT = Float64Col( pos=112 )
    p_CT_CGT = Float64Col( pos=113 )
    p_CT_ACGT = Float64Col( pos=114 )
    p_GG_AA = Float64Col( pos=115 )
    p_GG_AC = Float64Col( pos=116 )
    p_GG_AG = Float64Col( pos=117 )
    p_GG_AT = Float64Col( pos=118 )
    p_GG_CC = Float64Col( pos=119 )
    p_GG_CG = Float64Col( pos=120 )
    p_GG_CT = Float64Col( pos=121 )
    p_GG_GG = Float64Col( pos=122 )
    p_GG_GT = Float64Col( pos=123 )
    p_GG_TT = Float64Col( pos=124 )
    p_GG_ACG = Float64Col( pos=125 )
    p_GG_ACT = Float64Col( pos=126 )
    p_GG_AGT = Float64Col( pos=127 )
    p_GG_CGT = Float64Col( pos=128 )
    p_GG_ACGT = Float64Col( pos=129 )
    p_GT_AA = Float64Col( pos=130 )
    p_GT_AC = Float64Col( pos=131 )
    p_GT_AG = Float64Col( pos=132 )
    p_GT_AT = Float64Col( pos=133 )
    p_GT_CC = Float64Col( pos=134 )
    p_GT_CG = Float64Col( pos=135 )
    p_GT_CT = Float64Col( pos=136 )
    p_GT_GG = Float64Col( pos=137 )
    p_GT_GT = Float64Col( pos=138 )
    p_GT_TT = Float64Col( pos=139 )
    p_GT_ACG = Float64Col( pos=140 )
    p_GT_ACT = Float64Col( pos=141 )
    p_GT_AGT = Float64Col( pos=142 )
    p_GT_CGT = Float64Col( pos=143 )
    p_GT_ACGT = Float64Col( pos=144 )
    p_TT_AA = Float64Col( pos=145 )
    p_TT_AC = Float64Col( pos=146 )
    p_TT_AG = Float64Col( pos=147 )
    p_TT_AT = Float64Col( pos=148 )
    p_TT_CC = Float64Col( pos=149 )
    p_TT_CG = Float64Col( pos=150 )
    p_TT_CT = Float64Col( pos=151 )
    p_TT_GG = Float64Col( pos=152 )
    p_TT_GT = Float64Col( pos=153 )
    p_TT_TT = Float64Col( pos=154 )
    p_TT_ACG = Float64Col( pos=155 )
    p_TT_ACT = Float64Col( pos=156 )
    p_TT_AGT = Float64Col( pos=157 )
    p_TT_CGT = Float64Col( pos=158 )
    p_TT_ACGT = Float64Col( pos=159 )
    p_ACG_AA = Float64Col( pos=160 )
    p_ACG_AC = Float64Col( pos=161 )
    p_ACG_AG = Float64Col( pos=162 )
    p_ACG_AT = Float64Col( pos=163 )
    p_ACG_CC = Float64Col( pos=164 )
    p_ACG_CG = Float64Col( pos=165 )
    p_ACG_CT = Float64Col( pos=166 )
    p_ACG_GG = Float64Col( pos=167 )
    p_ACG_GT = Float64Col( pos=168 )
    p_ACG_TT = Float64Col( pos=169 )
    p_ACG_ACG = Float64Col( pos=170 )
    p_ACG_ACT = Float64Col( pos=171 )
    p_ACG_AGT = Float64Col( pos=172 )
    p_ACG_CGT = Float64Col( pos=173 )
    p_ACG_ACGT = Float64Col( pos=174 )
    p_ACT_AA = Float64Col( pos=175 )
    p_ACT_AC = Float64Col( pos=176 )
    p_ACT_AG = Float64Col( pos=177 )
    p_ACT_AT = Float64Col( pos=178 )
    p_ACT_CC = Float64Col( pos=179 )
    p_ACT_CG = Float64Col( pos=180 )
    p_ACT_CT = Float64Col( pos=181 )
    p_ACT_GG = Float64Col( pos=182 )
    p_ACT_GT = Float64Col( pos=183 )
    p_ACT_TT = Float64Col( pos=184 )
    p_ACT_ACG = Float64Col( pos=185 )
    p_ACT_ACT = Float64Col( pos=186 )
    p_ACT_AGT = Float64Col( pos=187 )
    p_ACT_CGT = Float64Col( pos=188 )
    p_ACT_ACGT = Float64Col( pos=189 )
    p_AGT_AA = Float64Col( pos=190 )
    p_AGT_AC = Float64Col( pos=191 )
    p_AGT_AG = Float64Col( pos=192 )
    p_AGT_AT = Float64Col( pos=193 )
    p_AGT_CC = Float64Col( pos=194 )
    p_AGT_CG = Float64Col( pos=195 )
    p_AGT_CT = Float64Col( pos=196 )
    p_AGT_GG = Float64Col( pos=197 )
    p_AGT_GT = Float64Col( pos=198 )
    p_AGT_TT = Float64Col( pos=199 )
    p_AGT_ACG = Float64Col( pos=200 )
    p_AGT_ACT = Float64Col( pos=201 )
    p_AGT_AGT = Float64Col( pos=202 )
    p_AGT_CGT = Float64Col( pos=203 )
    p_AGT_ACGT = Float64Col( pos=204 )
    p_CGT_AA = Float64Col( pos=205 )
    p_CGT_AC = Float64Col( pos=206 )
    p_CGT_AG = Float64Col( pos=207 )
    p_CGT_AT = Float64Col( pos=208 )
    p_CGT_CC = Float64Col( pos=209 )
    p_CGT_CG = Float64Col( pos=210 )
    p_CGT_CT = Float64Col( pos=211 )
    p_CGT_GG = Float64Col( pos=212 )
    p_CGT_GT = Float64Col( pos=213 )
    p_CGT_TT = Float64Col( pos=214 )
    p_CGT_ACG = Float64Col( pos=215 )
    p_CGT_ACT = Float64Col( pos=216 )
    p_CGT_AGT = Float64Col( pos=217 )
    p_CGT_CGT = Float64Col( pos=218 )
    p_CGT_ACGT = Float64Col( pos=219 )
    p_ACGT_AA = Float64Col( pos=220 )
    p_ACGT_AC = Float64Col( pos=221 )
    p_ACGT_AG = Float64Col( pos=222 )
    p_ACGT_AT = Float64Col( pos=223 )
    p_ACGT_CC = Float64Col( pos=224 )
    p_ACGT_CG = Float64Col( pos=225 )
    p_ACGT_CT = Float64Col( pos=226 )
    p_ACGT_GG = Float64Col( pos=227 )
    p_ACGT_GT = Float64Col( pos=228 )
    p_ACGT_TT = Float64Col( pos=229 )
    p_ACGT_ACG = Float64Col( pos=230 )
    p_ACGT_ACT = Float64Col( pos=231 )
    p_ACGT_AGT = Float64Col( pos=232 )
    p_ACGT_CGT = Float64Col( pos=233 )
    p_ACGT_ACGT = Float64Col( pos=234 )
