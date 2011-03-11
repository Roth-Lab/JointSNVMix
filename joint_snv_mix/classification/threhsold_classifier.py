'''
Created on 2011-02-10

@author: Andrew Roth
'''
import csv
import numpy as np

from joint_snv_mix import constants
from joint_snv_mix.classification.data import JointData
from joint_snv_mix.file_formats.jcnt import JointCountsReader

def run_threshold(args):
    runner = ThresholdRunner(args)
    
    runner.run(args)

#=======================================================================================================================
# Runner
#=======================================================================================================================
class ThresholdRunner(object):
    def __init__(self, args):
        self.model = ThresholdModel(args)
        
        self.data_class = JointData
        
        self.classes = ('Reference', 'Somatic', 'Somatic',
                        'LOH', 'Germline', 'LOH',
                        'Unknown', 'Unknown', 'Germline') 
    
    def run(self, args):
        self.reader = JointCountsReader(args.jcnt_file_name)
        self.writer = csv.writer(open(args.tsv_file_name, 'w'), delimiter='\t')
        
        chr_list = self.reader.get_chr_list()
        
        for chr_name in sorted(chr_list):
            self._classify_chromosome(chr_name)
                            
        self.reader.close()
        
    def _classify_chromosome(self, chr_name):
        counts = self.reader.get_counts(chr_name)
        jcnt_table = self.reader.get_table(chr_name)
        
        end = self.reader.get_chr_size(chr_name)

        n = int(1e5)
        start = 0
        stop = min(n, end)
        
        while start < end:
            sub_counts = counts[start:stop]
            sub_rows = jcnt_table[start:stop]
                              
            data = self.data_class(sub_counts)            
                
            labels = self.model.classify(data)
            
            self._write_rows(chr_name, sub_rows, labels)
        
            start = stop
            stop = min(stop + n, end)

    def _write_rows(self, chr_name, rows, labels):
        for i, row in enumerate(rows):
            out_row = [chr_name]
            out_row.extend(row)
            
            label = int(labels[i])
            
            class_name = self.classes[label]
                        
            out_row.append(class_name)
            
            if class_name == 'Somatic':
                print out_row
            
            self.writer.writerow(out_row)

#=======================================================================================================================
# Model
#=======================================================================================================================
class ThresholdModel(object):
    def __init__(self, args):        
        self.threshold = {}
        self.threshold['normal'] = args.normal_threshold
        self.threshold['tumour'] = args.tumour_threshold

        self.min_var_depth = args.min_var_depth
    
    def classify(self, data):
        genotypes = self._call_genotypes(data)
        
        joint_genotypes = self._call_joint_genotypes(data, genotypes)
        
        return joint_genotypes
    
    def _call_genotypes(self, data):
        genotypes = {}

        n = data.a['normal'].size
        
        for genome in constants.genomes:
            a = data.a[genome]
            b = data.b[genome]
            d = np.asanyarray(a + b, dtype=np.float)
            
            freq = b / d
            t = self.threshold[genome]
            
            ref = (freq < t)
            
            het = np.logical_and(t <= freq, freq <= (1 - t))
            
            hom = ((1 - t) < freq)
            
            genotypes[genome] = -1 * np.ones((n,))
            
            genotypes[genome][ref] = 0
            genotypes[genome][het] = 1
            genotypes[genome][hom] = 2
                    
        return genotypes
                    
    def _call_joint_genotypes(self, data, genotypes):
        normal = genotypes['normal']
        tumour = genotypes['tumour']
        
        joint_genotypes = 3 * normal + tumour
        
        b_T = data.b['tumour']
        
        # Set low coverage sites to reference.
        joint_genotypes[b_T < self.min_var_depth] = 0

        return joint_genotypes
