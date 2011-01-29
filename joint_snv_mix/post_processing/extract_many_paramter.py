import glob
import os
import sys

from extract_paramters import main as extract_params

dir = sys.argv[1]

glob_str = os.path.join( dir, '*.jsm' )

file_list = glob.glob( glob_str )

for file_name in file_list:
    print "#" * 100
    print file_name
    print "#" * 100
    
    extract_params( file_name )
    
    print
    print
    
