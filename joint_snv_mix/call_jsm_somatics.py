import bisect
import csv

from joint_snv_mix.file_formats.jsm import JointSnvMixReader

fields = [
          'chrom',
          'position',
          'ref_base',
          'normal_var_base',
          'tumour_var_base',
          'normal_counts_a',
          'normal_counts_b',
          'tumour_counts_a',
          'tumour_counts_b',
          'somatic_prob'          
          ]

genotypes = ['aa', 'ab', 'bb']

for g_1 in genotypes:
    for g_2 in genotypes:        
        prob_field = "_".join(('p', g_1, g_2))
        fields.append(prob_field)

def jsm_to_tsv(args):
    '''
    Call somatic mutations form jsm file and output to human readable tab separated file in sorted order descending by
    somatic probability.
    '''
    jsm_file_name = args.jsm_file_name
    out_file_name = args.out_file_name    

    reader = JointSnvMixReader(jsm_file_name)
    writer = csv.DictWriter(open(out_file_name, 'w'), fields, delimiter='\t')
    writer.writeheader()
    
    table_list = reader.get_table_list()
    
    for chr_name in sorted(table_list):
        table = reader.get_table(chr_name)
        
        for row in table:
            row = format_rows(row, chr_name)            
            writer.writerow(row)
            
    reader.close()

def call_somatics_from_jsm(args):
    '''
    Call somatic mutations form jsm file and output to human readable tab separated file in sorted order descending by
    somatic probability.
    '''
    jsm_file_name = args.jsm_file_name
    out_file_name = args.out_file_name
    
    if args.auto:
        rows = load_auto_threshold_somatics(jsm_file_name)
    else:
        rows = load_manual_threshold_somatics(jsm_file_name, args.threshold)

    rows.reverse()

    writer = csv.DictWriter(open(out_file_name, 'w'), fields, delimiter='\t')

    writer.writeheader()
    writer.writerows(rows)
    
def load_manual_threshold_somatics(jsm_file_name, threshold):
    '''
    Load a list of rows containing somatics based on pre-specified probability threshold.
    '''
    reader = JointSnvMixReader(jsm_file_name)

    chr_list = reader.get_table_list()

    rows = []
    scores = []

    for chr_name in sorted(chr_list):
        print chr_name

        chr_rows = reader.get_table(chr_name)

        for row in chr_rows:
            score = row['p_aa_ab'] + row['p_aa_bb']            
            
            if score >= threshold:
                row = format_rows(row, chr_name)
                
                insert_position = bisect.bisect(scores, score)
                
                scores.insert(insert_position, score)
                rows.insert(insert_position, row)
                
    reader.close()

    return rows


def load_auto_threshold_somatics(jsm_file_name):
    '''
    Load a list of rows containing somatics based on automatically determined probability threshold. Threshold is
    determined based on inflection point method.
    '''
    n = int(1e5)
    threshold = 1e-6

    reader = JointSnvMixReader(jsm_file_name)

    chr_list = reader.get_table_list()

    scores = []
    rows = []

    for chr_name in sorted(chr_list):
        print chr_name

        chr_rows = reader.get_table(chr_name)

        for row in chr_rows:
            score = row['p_aa_ab'] + row['p_aa_bb']

            insert_position = bisect.bisect(scores, score)

            if insert_position > 0 or len(scores) == 0:
                scores.insert(insert_position, score)
                
                row = format_rows(row, chr_name)
                
                rows.insert(insert_position, row)
            
                if scores[0] <= threshold or len(scores) > n:
                    scores.pop(0)
                    rows.pop(0)

    reader.close()
    
    max_diff = 0
    index = 0
    
    for i in range(len(scores) - 1):
        diff = scores[i + 1] - scores[i]
        
        if diff > max_diff:
            max_diff = diff
            index = i
            
    rows = rows[index:]

    return rows

def format_rows(row, chrom):
    row = row.fetch_all_fields()
    row = dict(zip(row.dtype.names, row.real))
    row['somatic_prob'] = row['p_aa_ab'] + row['p_aa_bb']
    row['chrom'] = chrom
    
    row['normal_var_base'] = row['normal_base']
    row['tumour_var_base'] = row['tumour_base']
    
    del row['normal_base']
    del row['tumour_base']
    
    return row

if __name__ == "__main__":
    import sys
    
    jsm_file_name = sys.argv[1]
    
    out_file_name = sys.argv[2]
    
    call_somatics_from_jsm(jsm_file_name, out_file_name)
