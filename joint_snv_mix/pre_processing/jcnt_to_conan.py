'''
Created on 2011-02-10

@author: Andrew Roth
'''
import csv

import numpy as np

from joint_snv_mix.file_formats.jcnt import JointCountsReader
from joint_snv_mix.file_formats.cncnt import ConanCountsFile
from argparse import Namespace

def jcnt_to_cncnt(args):
    reader = JointCountsReader(args.jcnt_file_name)
    
    chr_list = reader.get_chr_list()
    
    cncnt_file = ConanCountsFile(args.cncnt_file_name, 'w')
    
    segment_reader = csv.reader(open(args.segment_file_name), delimiter='\t')

    old_chr_name = None

    for row in segment_reader:
        chr_name = row[0]
        
        if chr_name == '23':
            chr_name = 'X'
        if chr_name == '24':
            chr_name = 'Y'
        
        start = int(row[1])
        stop = int(row[2])
        cn_status = row[3]
        
        if cn_status == '7':
            cn_status = '1'
        elif cn_status == '8':
            cn_status = '2'
        elif cn_status == '9':
            cn_status = '4'
        elif cn_status == '10':
            cn_status = '5'
        elif cn_status == '11':
            cn_status = '6'
        
        if chr_name not in chr_list:
            print chr_name
            continue
        
        # Load new table if chromosome has changed.
        if chr_name != old_chr_name:
            rows = reader.get_table(chr_name)
        
        old_chr_name = chr_name        
        
        segment_indices = np.logical_and(rows['position'] >= start, rows['position'] <= stop)
        
        segment_rows = rows[segment_indices]
        
        if len(segment_rows) == 0:
            continue
        
        cncnt_file.add_rows(cn_status, chr_name, segment_rows)
    
    reader.close()
    cncnt_file.close()
        
if __name__ == "__main__":
    import sys
    
    args = Namespace()
    
    args.jcnt_file_name = sys.argv[1]
    args.cncnt_file_name = sys.argv[2]    
    args.segment_file_name = sys.argv[3]
    
    
    jcnt_to_cncnt(args)
