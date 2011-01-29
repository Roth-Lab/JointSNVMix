'''
Created on 2010-07-29

Classes for generating fake data for unit testing.

@author: Andrew Roth
'''
import random

import numpy as np

from pyleup.pileup.binary_pileup import ReadsColumn

chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
            '11', '12', '13', '14', '15', '16', '17', '18', '19',
            '20', '21', '22', 'X', 'Y']

alphabet = ['A', 'C', 'G', 'T']

unknown_alphabet = ['N', '*']

class ReadsColumnGenerator:
    def create_read( self, length, unknown_freq = 0.01 ):
        chr_name = 'chr' + random.choice( chr_list )

        position = random.randint( 0, int( 1e9 ) )

        ref_base = random.choice( alphabet )

        bases = []

        for i in range( length ):
            bases.append( self._get_base( unknown_freq ) )

        reads = np.array( bases, dtype = [( 'base', 'a1' ), ( 'base_qual', np.uint8 ), ( 'map_qual', np.uint8 )] )

        reads = reads.view( np.recarray )

        column = ReadsColumn( chr_name, position, ref_base, reads )

        expected_values = {}
        expected_values['chr_name'] = chr_name
        expected_values['position'] = position
        expected_values['ref_base'] = ref_base
        expected_values['bases'] = bases

        return expected_values, column



    def _get_base( self, unknown_freq ):
        if random.random() <= unknown_freq:
            base = random.choice( unknown_alphabet )
        else:
            base = random.choice( alphabet )

        base_qual = random.randint( 0, 93 )

        map_qual = random.randint( 0, 93 )

        return ( base, base_qual, map_qual )

    def get_counts( self, bases, min_base_qual = 0, min_map_qual = 0 ):
        counts = {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}

        for base_read in bases:
            base = base_read[0]
            base_qual = base_read[1]
            map_qual = base_read[2]

            if base in counts and base_qual >= min_base_qual and map_qual >= min_map_qual:
                counts[base] += 1
        
        return counts

