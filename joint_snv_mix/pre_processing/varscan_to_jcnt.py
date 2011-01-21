import csv
import sys

from joint_snv_mix.file_formats.jcnt import JointCountsFile

def write_rows( jcnt_file, rows ):
    for chr_name, chr_rows in rows.items():
        jcnt_file.add_rows( chr_name, chr_rows )

#=======================================================================================================================
# Top Level Code
#=======================================================================================================================
varscan_file_name = sys.argv[1]
jcnt_file_name = sys.argv[2]

varscan_reader = csv.DictReader( open( varscan_file_name ), delimiter='\t' )
jcnt_file = JointCountsFile( jcnt_file_name, 'w' )

nucs = ['A', 'C', 'G', 'T']
    
rows = {}
i = 0
    
for row in varscan_reader:
    # Skip rows without a variant in the tumour or normal or with indel or >= 3 bases.
    if row['var'] not in nucs:
        continue
    
    chr_name = row['chrom']
    
    if chr_name not in rows:
        rows[chr_name] = []
    
    chr_coord = int( row['position'] )
    
    normal_ref_counts = int( row['normal_reads1'] )
    tumour_ref_counts = int( row['tumor_reads1'] )
    
    normal_non_ref_counts = int( row['normal_reads2'] )
    tumour_non_ref_counts = int( row['tumor_reads2'] )
    
    normal_counts = ( normal_ref_counts, normal_non_ref_counts )
    tumour_counts = ( tumour_ref_counts, tumour_non_ref_counts )
    
    ref_base = row['ref']
    
    normal_non_ref_base = row['var']
    tumour_non_ref_base = row['var']
    
    jcnt_entry = [ chr_coord, ref_base, normal_non_ref_base, tumour_non_ref_base ]
    
    jcnt_entry.extend( normal_counts )
    jcnt_entry.extend( tumour_counts )
    
    rows[chr_name].append( jcnt_entry )
    
    print "\t".join( [str( x ) for x in jcnt_entry] )
    
    i += 1
    
    if i >= 1e4:
        print chr_name, chr_coord
        
        write_rows( jcnt_file, rows )
        
        rows = {}
        i = 0
    
# Last call to write remaining rows.
write_rows( jcnt_file, rows )
    
jcnt_file.close()
