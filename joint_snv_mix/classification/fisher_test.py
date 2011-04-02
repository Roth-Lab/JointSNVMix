'''
Created on 2011-02-09

@author: Andrew Roth
'''
import csv

import numpy as np

from fisher import pvalue_npy

from joint_snv_mix import constants

from joint_snv_mix.file_formats.jcnt import JointCountsReader

from joint_snv_mix.classification.utils.data import JointData

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
    def __init__(self, p_value_threshold=0.05, base_line_error=0.001,
                 min_var_freq=0.1, min_hom_freq=0.9, min_var_depth=4):
        
        self.p_value_threshold = p_value_threshold
        
        self.base_line_error = base_line_error
        
        self.min_var_freq = min_var_freq
        self.min_hom_freq = min_hom_freq
        self.min_var_depth = min_var_depth
        
    def classify(self, data):
        genotypes = self._call_genotypes(data)
        
        joint_genotypes = self._call_joint_genotypes(data, genotypes)
        
        return joint_genotypes
    
    def _get_right_tail_significance(self, a, b, expected_freq):
        '''
        Compute p-value of a and b given null hypothesis b is due to random error with expected_freq as error rate.
        '''
        d = a + b
        
        expected_b = np.around(d * expected_freq)
        
        expected_a = d - expected_b
        
        # Downcast to uint to work with fisher exact test function.
        a = np.asarray(a, dtype=np.uint)
        b = np.asarray(b, dtype=np.uint)        
        expected_a = np.asarray(expected_a, dtype=np.uint)
        expected_b = np.asarray(expected_b, dtype=np.uint)
        
        left_tail, right_tail, two_tail = self._get_pvalues(expected_a, expected_b, a, b)
        
        return right_tail
    
    def _get_two_tail_significance(self, a_N, b_N, a_T, b_T):
        left_tail, right_tail, two_tail = self._get_pvalues(a_N, b_N, a_T, b_T)
        
        return two_tail
        
    def _get_pvalues(self, a_1, b_1, a_2, b_2):
        a_1 = np.asarray(a_1, dtype=np.uint)
        b_1 = np.asarray(b_1, dtype=np.uint)        
        a_2 = np.asarray(a_2, dtype=np.uint)
        b_2 = np.asarray(b_2, dtype=np.uint)
        
        left_tail, right_tail, two_tail = pvalue_npy(
                                                     a_1,
                                                     b_1,
                                                     a_2,
                                                     b_2
                                                     )
        
        return left_tail, right_tail, two_tail

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
                        
            p_values = self._get_right_tail_significance(a, b, self.base_line_error)
                        
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

        joint_genotypes = 3 * normal + tumour
        joint_genotypes = np.asarray(joint_genotypes, dtype=np.int)
        
        return joint_genotypes

class JointFisherModel(FisherModel):                        
    def _call_joint_genotypes(self, data, genotypes):
        normal = genotypes['normal']
        tumour = genotypes['tumour']
        
        joint_genotypes = 3 * normal + tumour
        joint_genotypes = np.asarray(joint_genotypes, dtype=np.int)
        
        a_N = data.a['normal']
        b_N = data.b['normal']
        a_T = data.a['tumour']
        b_T = data.b['tumour']
        
        # Check if the counts in the tumour are significantly different from normal using two tailed fisher.
        p_values = self._get_two_tail_significance(a_N, b_N, a_T, b_T)
        
        # Find places where the called genotypes do not agree in tumour normal and check significance.
        non_match = (normal != tumour)
        significant_p_values = (p_values <= self.p_value_threshold)                
        non_significant_p_values = np.logical_not(significant_p_values)        
        non_significant_non_match = np.logical_and(non_match, non_significant_p_values)
        
        # Change all calls where genotypes in tumour/normal differ based on independent calls, but whose counts do not
        # differ significantly to the genotype (AB,AB) i.e. het germline.
        joint_genotypes[non_significant_non_match] = 4
        
        return joint_genotypes
