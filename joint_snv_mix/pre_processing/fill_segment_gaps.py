'''
Created on 2011-02-11

@author: Andrew
'''
import csv

def main(args):
    unknown_class_id = 'unknown'
    
    reader = csv.reader(open(args.seg_file), delimiter='\t')
    writer = csv.writer(open(args.out_file, 'w'), delimiter='\t')
    
    row = reader.next()
        
    chr_name = row[0]
    start = int(row[1])
    stop = int(row[2])
    
    if start != 0:
        out_row = [chr_name, 0, start - 1, unknown_class_id]
        writer.writerow(out_row)
        
    writer.writerow(row)
    
    for row in reader:
        old_chr_name = chr_name
        old_stop = stop

        chr_name = row[0]
        start = int(row[1])
        stop = int(row[2])
        
        if chr_name != old_chr_name:
            # Close old chromosome out.
            out_row = [old_chr_name, old_stop + 1, int(1e10), unknown_class_id]            
            writer.writerow(out_row)
            
            # Open new chromosome
            if start != 0:
                out_row = [chr_name, 0, start - 1, unknown_class_id]
        
                writer.writerow(out_row)
                
                old_stop = start - 1
        
        # Fill gap in middle of chromosome        
        if old_stop != start - 1:
            out_row = [chr_name, old_stop + 1, start - 1, unknown_class_id]
            writer.writerow(out_row)
        
        writer.writerow(row)
    
    # Close last chromosome out.
    out_row = [chr_name, stop + 1, int(1e10), unknown_class_id]            
    writer.writerow(out_row)
        

if __name__ == "__main__":
    import sys
    from argparse import Namespace
    
    args = Namespace()
    
    args.seg_file = sys.argv[1]
    args.out_file = sys.argv[2]
    
    main(args)
