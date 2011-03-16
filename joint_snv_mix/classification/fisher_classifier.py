'''
Created on 2011-02-09

@author: Andrew Roth
'''
import csv

import numpy as np

from fisher import pvalue_npy

from joint_snv_mix import constants
from joint_snv_mix.classification.data import JointData
from joint_snv_mix.file_formats.jcnt import JointCountsReader


def run_fisher(args):            
    if args.model == "joint":
        runner = JointFisherRunner(args)
        
    elif args.model == "independent":
        runner = IndependentFisherRunner(args)
        
    runner.run(args)
 
#=======================================================================================================================
# Runner Code
#=======================================================================================================================
class FisherRunner(object):
    def __init__(self):
        self.data_class = JointData
        
        self.classes = ('Reference', 'Germline', 'Somatic', 'LOH', 'Unknown')
    
    def run(self, args):
        self.reader = JointCountsReader(args.jcnt_file_name)
        self.writer = csv.writer(open(args.tsv_file_name, 'w'), delimiter='\t')
        
        chr_list = self.reader.get_table_list()
        
        for chr_name in sorted(chr_list):
            self._classify_chromosome(chr_name)
                            
        self.reader.close()
        
    def _classify_chromosome(self, chr_name):
        counts = self.reader.get_counts(chr_name)
        jcnt_table = self.reader.get_table(chr_name)
        
        end = self.reader.get_number_of_table_rows(chr_name)

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

class IndependentFisherRunner(FisherRunner):
    def __init__(self, args):
        self.model = IndependentFisherModel(args)
        
        FisherRunner.__init__(self)
        
class JointFisherRunner(FisherRunner):
    def __init__(self, args):
        self.model = JointFisherModel(args)
        
        FisherRunner.__init__(self)

#=======================================================================================================================
# Models
#=======================================================================================================================
class FisherModel(object):
    def __init__(self, args):
        
        self.p_value_threshold = args.p_value_threshold
        
        self.base_line_error = args.base_line_error
        
        self.min_var_freq = args.min_var_freq
        self.min_hom_freq = args.min_hom_freq
        self.min_var_depth = args.min_var_depth
        
    def classify(self, data):
        genotypes = self._call_genotypes(data)
        
        joint_genotypes = self._call_joint_genotypes(data, genotypes)
        
        return joint_genotypes
    
    def _get_significance(self, a, b, expected_freq):
        '''
        Compute the p-value where the null hypothesis is that b was obtained due to error.
        
        a, b : numpy array of counts
        base_line_error : expected error rate.
        '''
        d = a + b
        
        expected_b = np.around(d * expected_freq)
        
        expected_a = d - expected_b
        
        # Downcast to uint to work with fisher exact test function.
        a = np.asarray(a, dtype=np.uint)
        b = np.asarray(b, dtype=np.uint)
        expected_a = np.asarray(expected_a, dtype=np.uint)
        expected_b = np.asarray(expected_b, dtype=np.uint)
        
        left_tail, right_tail, two_tail = pvalue_npy(expected_a,
                                                      expected_b,
                                                      a,
                                                      b)
        
        return right_tail

    def _call_joint_genotypes(self, data, genotypes):
        '''
        Return 0 = ref, 1 = germline, 2 = somatic, 3 = loh, 4 = unkown
        '''
        raise NotImplemented
    
    def _call_genotypes(self, data):
        genotypes = {}

        n = data.a['normal'].size
        
        for genome in constants.genomes:
            a = data.a[genome]
            b = data.b[genome]
            d = np.asanyarray(a + b, dtype=np.float)
                        
            p_values = self._get_significance(a, b, self.base_line_error)
                        
            below_threshold = (p_values <= self.p_value_threshold)
            
            var_freq = b / d            

            # Find variants below p-value threshold and above minimum variant frequency.
            above_min_var_freq = (var_freq >= self.min_var_freq)
            
            var_indices = np.logical_and(below_threshold, above_min_var_freq)
            
            # Call homozygous if above minimum threshold.
            above_min_hom_freq = (var_freq >= self.min_hom_freq)                        
            hom_var_indices = np.logical_and(var_indices, above_min_hom_freq)
            
            not_above_min_hom_freq = np.logical_not(above_min_hom_freq)
            het_var_indices = np.logical_and(var_indices, not_above_min_hom_freq)
            
            # Apply min variant depth filter
            above_min_var_depth = (b >= self.min_var_depth)
            
            het_var_indices = np.logical_and(het_var_indices, above_min_var_depth)
            hom_var_indices = np.logical_and(hom_var_indices, above_min_var_depth)
            
            # Assign genotypes assuming reference as default
            genotypes[genome] = np.zeros((n,))
            
            genotypes[genome][het_var_indices] = 1
            genotypes[genome][hom_var_indices] = 2
                    
        return genotypes

class IndependentFisherModel(FisherModel):
    def _call_joint_genotypes(self, data, genotypes):
        normal = genotypes['normal']
        tumour = genotypes['tumour']
                
        normal_aa = (normal == 0)
        normal_ab = (normal == 1)
        normal_bb = (normal == 2)
            
        normal_var = np.logical_or(normal_ab, normal_bb)
        
        tumour_aa = (tumour == 0)
        tumour_ab = (tumour == 1)
        tumour_bb = (tumour == 2)
        
        tumour_var = np.logical_or(tumour_ab, tumour_bb)
        tumour_hom = np.logical_and(tumour_aa, tumour_bb)
        
        reference = np.logical_and(normal_aa, tumour_aa)
        germline = np.logical_and(normal_var, tumour_var)
        somatic = np.logical_and(normal_aa, tumour_var)
        loh = np.logical_and(normal_ab, tumour_hom)
        
        
        n = normal_aa.size
        joint_genotypes = 4 * np.ones((n,))
        
        joint_genotypes[reference] = 0
        joint_genotypes[germline] = 1
        joint_genotypes[somatic] = 2
        joint_genotypes[loh] = 3
        
        return joint_genotypes

class JointFisherModel(FisherModel):                        
    def _call_joint_genotypes(self, data, genotypes):
        normal = genotypes['normal']
        tumour = genotypes['tumour']
        
        normal_aa = (normal == 0)
        normal_ab = (normal == 1)
        normal_bb = (normal == 2)
            
        normal_var = np.logical_or(normal_ab, normal_bb)
        
        tumour_aa = (tumour == 0)
        tumour_ab = (tumour == 1)
        tumour_bb = (tumour == 2)
        
        tumour_var = np.logical_or(tumour_ab, tumour_bb)
                
        reference = np.logical_and(normal_aa, tumour_aa)
        germline = np.logical_and(normal_var, tumour_var)
        
        a_N = data.a['normal']
        b_N = data.b['normal']
        d_N = np.asanyarray(a_N + b_N, dtype=np.float)      
        normal_freq = b_N / d_N
        
        a_T = data.a['tumour']
        b_T = data.b['tumour']
        
        p_values = self._get_significance(a_T, b_T, normal_freq)
        
        non_match = (normal != tumour)
        significant_p_values = (p_values <= self.p_value_threshold)        
        significant_non_match = np.logical_and(non_match, significant_p_values)
        
        somatic = np.logical_and(normal_aa, significant_non_match)
        loh = np.logical_and(normal_ab, significant_non_match)
        
        uknown = np.logical_and(normal_bb, significant_non_match)
        
        non_significant_p_values = np.logical_not(significant_p_values)
        non_significant_non_match = np.logical_and(non_match, non_significant_p_values)        
        germline = np.logical_or(germline, non_significant_non_match)
               
        n = a_N.size
        joint_genotypes = -1 * np.ones((n,))
        
        joint_genotypes[reference] = 0
        joint_genotypes[germline] = 1
        joint_genotypes[somatic] = 2
        joint_genotypes[loh] = 3
        joint_genotypes[uknown] = 4
        
        return joint_genotypes
