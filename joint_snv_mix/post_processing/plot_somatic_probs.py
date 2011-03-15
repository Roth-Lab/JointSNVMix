import matplotlib
matplotlib.use( 'agg' )

import matplotlib.pyplot as plot

from joint_snv_mix.file_formats.jsm import JointSnvMixReader

import bisect

excluded_chrom = ['Y', 'MT']

def main( jsm_file_name, out_file_name ):
    somatics = load_somatics( jsm_file_name )

    fig = plot.figure()

    ax = fig.add_subplot( 1, 1, 1 )

    ax.set_ylim( -0.05, 1.05 )

    ax.plot( somatics )

    plot.savefig( out_file_name )


def load_somatics( jsm_file_name ):
    n = int( 1e5 )
    threshold = 1e-6

    reader = JointSnvMixReader( jsm_file_name )

    chr_list = reader.get_chr_list()

    scores = []

    for chr_name in sorted( chr_list ):
        if chr_name in excluded_chrom:
            continue
        
        print chr_name

        chr_rows = reader.get_table( chr_name )

        for row in chr_rows:
            score = row['p_aa_ab'] + row['p_aa_bb']

            insert_position = bisect.bisect( scores, score )

            if insert_position > 0 or len( scores ) == 0:
                scores.insert( insert_position, score )
            
                if scores[0] <= threshold or len( scores ) > n:
                    scores.pop( 0 )

    reader.close()
    
    max_diff = 0
    index = 0
    
    for i in range( len( scores ) - 1 ):
        diff = scores[i + 1] - scores[i]
        
        if diff > max_diff:
            max_diff = diff
            index = i
            
    scores = scores[index:]

    return scores

if __name__ == "__main__":
    import sys
    
    jsm_file_name = sys.argv[1]
    
    out_file_name = sys.argv[2]
    
    main( jsm_file_name, out_file_name )
