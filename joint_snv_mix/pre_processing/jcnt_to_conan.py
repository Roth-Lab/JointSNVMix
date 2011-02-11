'''
Created on 2011-02-10

@author: Andrew Roth
'''
import csv

import numpy as np

from joint_snv_mix.file_formats.jcnt import JointCountsReader
from joint_snv_mix.file_formats.cncnt import ConanCountsFile
from argparse import Namespace


def main( args ):
    reader = JointCountsReader( args.jcnt_file_name )
    
    chr_list = reader.get_chr_list()
    
    cncnt_file = ConanCountsFile( args.cncnt_file_name, 'w' )
    
    segment_reader = csv.reader( open( args.segment_file_name ), delimiter='\t' )

    for row in segment_reader:
        print row
        
        chr_name = row[0]
        start = int( row[1] )
        stop = int( row[2] )
        cn_status = row[3]
        
        if chr_name not in chr_list:
            continue
        
        rows = reader.get_rows( chr_name )
        
        segment_indices = np.logical_and( rows['position'] >= start, rows['position'] <= stop )
        
        cncnt_file.add_rows( cn_status, chr_name, rows[segment_indices] )

    
    reader.close()
    cncnt_file.close()
        
if __name__ == "__main__":
    import sys
    
    args = Namespace()
    
    args.jcnt_file_name = sys.argv[1]    
    args.segment_file_name = sys.argv[2]
    args.cncnt_file_name = sys.argv[3]
    
    main( args )
