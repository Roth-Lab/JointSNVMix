'''
Created on 2012-01-26

@author: Andrew Roth
'''
from collections import defaultdict

from joint_snv_mix.counter import JointBinaryCounter
from joint_snv_mix.positions_counter import PositionsCounter
from joint_snv_mix.post_processing.feature_extractor import FeatureExtractor
from joint_snv_mix.samtools import BamFile, FastaFile

import csv
import os
import tempfile

def main(ref_genome_file, case_list_file, ground_truth_file, out_file_name):
    ref_genome = FastaFile(ref_genome_file)
    
    training_files = get_training_files(case_list_file)
    
    features = defaultdict(dict)
    labels = {}
    
    for case in training_files:
        print case
        labels[case] = get_labels(case, ground_truth_file)
        
        training_files[case]['positions_file'] = get_positions_file(case, ground_truth_file) 
        
        features[case] = extract_features(ref_genome, training_files[case])
        
        os.unlink(training_files[case]['positions_file'])
    
    write_features(features, labels, out_file_name)

def get_labels(case, ground_truth_file):
    labels = {}
    
    print ground_truth_file
    
    reader = csv.DictReader(open(ground_truth_file), delimiter='\t')
    
    for row in reader:
        if row['case'] != case:
            continue
    
        pos = (row['chrom'], int(row['coord']))
        
        labels[pos] = int(row['label'])
    
    return labels

def get_training_files(case_list_file):
    training_files = defaultdict(dict)
    
    reader = csv.DictReader(open(case_list_file), delimiter='\t')
    
    for row in reader:
        case = row['case']

        training_files[case]['normal_file'] = row['normal_file']
        training_files[case]['tumour_file'] = row['tumour_file']
    
    return training_files

def get_positions_file(case, ground_truth_file):
    temp_file, file_name = tempfile.mkstemp()
    
    print temp_file, file_name
    
    os.close(temp_file)
        
    reader = csv.DictReader(open(ground_truth_file), delimiter='\t')
    writer = csv.writer(open(file_name,'w'), delimiter=' ')
    
    positions = []
    
    for row in reader:
        if case != row['case']:
            continue
        
        chrom = row['chrom']
        coord = int(row['coord'])
        
        pos = (chrom, coord)
        
        positions.append(pos)
    
    writer.writerows(sorted(positions))
    
    return file_name

def extract_features(ref_genome, case_info):
    features = {}
    
    normal_bam = BamFile(case_info['normal_file'])
    tumour_bam = BamFile(case_info['tumour_file'])
    
    extractor = FeatureExtractor(normal_bam, tumour_bam)
    
    counter = JointBinaryCounter(normal_bam, tumour_bam, ref_genome, qualities=1)
    
    positions_counter = PositionsCounter(case_info['positions_file'], counter)    
    
    for ref in positions_counter.refs:
        print ref
        ref_iter = positions_counter.get_ref_iterator(counter.get_ref_iterator(ref))
        
        for row in ref_iter:
            pos = (row.ref, row.position)
            features[pos] = extractor.get_features(row)
    
    return features

def write_features(features, labels, out_file_name):
    writer = csv.writer(open(out_file_name, 'w'), delimiter='\t')
    
    for case in features:
        for pos in features[case]:
            pos_id = "_".join((case, pos[0], str(pos[1])))
            
            out_row = [pos_id, ]
            
            out_row.extend(features[case][pos])
                        
            out_row.append(labels[case][pos])
            
            writer.writerow(out_row) 
            
if __name__ == "__main__":
    import sys

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
