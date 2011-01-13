'''
Created on 2010-12-09

@author: Andrew Roth
'''
import math
import random

import numpy as np

def subsample( reader, sample_size ):
    chr_list = reader.get_chr_list()
    
    sample = []
    
    nrows = reader.get_data_set_size()
    
    for chr_name in chr_list:
        chr_size = reader.get_chr_size( chr_name=chr_name )
        
        chr_sample_size = math.ceil( float( chr_size ) / nrows * sample_size )
        chr_sample_size = int( chr_sample_size )
        
        chr_sample_indices = random.sample( xrange( chr_size ), chr_sample_size )
        
        chr_counts = reader.get_counts( chr_name )
        
        chr_sample = chr_counts[chr_sample_indices]
        
        sample.append( chr_sample )
        
    sample = np.vstack( sample )
    
    return sample
