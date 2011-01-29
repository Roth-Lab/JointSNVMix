import matplotlib.pyplot as plot

from joint_snv_mix.file_formats.jsm import JointSnvMixReader

import bisect

excluded_chrom = ['Y', 'MT']

def main( jsm_file_name ):
    n = int( 1e5 )

    somatics = load_somatics( jsm_file_name, n )

    plot.plot( somatics )

    plot.savefig( 'somatic_probs.pdf' )


def load_somatics( jsm_file_name, n ):
    reader = JointSnvMixReader( jsm_file_name )

    chr_list = reader.get_chr_list()

    position_score = []
    scores = []

    for chr_name in sorted( chr_list ):
        if chr_name in excluded_chrom:
            continue
        chr_rows = reader.get_rows( chr_name )

        for row in chr_rows:
            position = int( row['position'] )
            score = float( row['p_aa_ab'] + row['p_aa_bb'] )

            insert_position = bisect.bisect( scores, score )

            if insert_position > 0 or len( scores ) == 0:
                position_score.insert( insert_position, ( chr_name, position, score ) )
                scores.insert( insert_position, score )

                if len( scores ) >= n:
                    scores.pop( 0 )
                    position_score.pop( 0 )

                print position, insert_position, position_score[0], position_score[-1]


    reader.close()

    return position_score
