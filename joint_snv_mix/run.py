'''
Created on 2011-08-11

@author: Andrew Roth
'''
import csv
import sys

from joint_snv_mix.counter import JointBinaryCounter
from joint_snv_mix.models.joint_snv_mix import JointSnvMixModel, JointSnvMixPriors, JointSnvMixParameters
from joint_snv_mix.samtools import BamFile, FastaFile

#=======================================================================================================================
# Functions for running classification.
#=======================================================================================================================
def classify(args):
    counter_factory = CounterFactory()
    model_factory = ModelFactory()

    counter = counter_factory.get_counter(args)
    model = model_factory.get_model(args) 
    
    classify_data_set(counter, model, args)

def classify_data_set(counter, classifier, args):
    writer = ResultsWriter(file_name=args.out_file)
    
    if args.chromosome is not None:
        positions_iter = counter.get_ref_iterator(args.chromosome)
    else:
        positions_iter = get_genome_iterator(counter)
    
    for row in positions_iter:
        if not args.print_all_positions and row.tumour_var_counts == 0:
            continue 
        
        probs = classifier.predict(row.data)
        
        writer.write_position(row, probs)

class ResultsWriter(object):
    info = [
            'chrom',
            'position',
            'ref_base',
            'var_base',
            'normal_counts_a',
            'normal_counts_b',
            'tumour_counts_a',
            'tumour_counts_b'
            ]
    
    probs = [
             'p_AA_AA',
             'p_AA_AB',
             'p_AA_BB',
             'p_AB_AA',
             'p_AB_AB',
             'p_AB_BB',
             'p_BB_AA',
             'p_BB_AB',
             'p_BB_BB'                            
             ]
    
    fields = info + probs
    
    def __init__(self, file_name=None):
        if file_name is None:
            fh = sys.stdout
        else:
            fh = open(file_name, 'w')
        
        self._writer = csv.DictWriter(fh, ResultsWriter.fields, delimiter='\t')
        
        self._writer.writeheader()
        
    def write_position(self, row, probs):
        out_row = {}
        out_row['chrom'] = row.ref
        out_row['position'] = row.position
        out_row['ref_base'] = row.ref_base
        out_row['var_base'] = row.var_base
        
        counts = row.counts
        
        out_row['normal_counts_a'] = counts[0]
        out_row['normal_counts_b'] = counts[1]
        out_row['tumour_counts_a'] = counts[2]
        out_row['tumour_counts_b'] = counts[3]
        
        probs = dict(zip(ResultsWriter.probs, probs))
        
        out_row = dict(out_row.items() + probs.items())
        
        self._writer.writerow(out_row)
         
#=======================================================================================================================
# Functions for training a model.
#=======================================================================================================================
def train(args):
    counter_factory = CounterFactory()
    model_factory = ModelFactory()
    
    counter = counter_factory.get_counter(args)
    model = model_factory.get_model(args)
    
    data_set = create_training_data_set(counter, args)
    
    model.fit(data_set, verbose=True)
    
    model.params.write_to_file(args.estimated_parameters_file)        

def create_training_data_set(counter, args):
    i = 0
    
    data_set = []
    
    if args.chromosome is not None:
        positions_iter = counter.get_ref_iterator(args.chromosome)
    else:
        positions_iter = get_genome_iterator(counter)
    
    for row in positions_iter:
        if args.min_normal_depth <= row.normal_depth <= args.max_normal_depth and \
           args.min_tumour_depth <= row.tumour_depth <= args.max_tumour_depth:
            i += 1
            
            if i % args.skip_size == 0:
                data_set.append(row.data)
    
    return data_set

#=======================================================================================================================
# Helper factory classes for setting up parsers and classifiers.
#=======================================================================================================================
def get_genome_iterator(counter):
    refs = sorted(counter.refs)
    
    for ref in refs:
        ref_iter = counter.get_ref_iterator(ref)
        
        for row in ref_iter:
            yield row

class CounterFactory(object):
    def get_counter(self, args):
        if args.model == 'jsm1':
            qualities = 0
            min_base_qual = args.min_base_qual
            min_map_qual = args.min_map_qual
            
        elif args.model == 'jsm2':
            qualities = 1
            min_base_qual = 0
            min_map_qual = 0
        
        normal_bam = BamFile(args.normal_file)
        tumour_bam = BamFile(args.tumour_file)
        ref_genome = FastaFile(args.reference_genome_file)
        
        
        counter = JointBinaryCounter(normal_bam,
                                     tumour_bam,
                                     ref_genome,
                                     min_base_qual,
                                     min_map_qual,
                                     qualities)
        
        return counter

class ModelFactory(object):
    def get_model(self, args):
        priors = self._get_joint_snv_mix_priors(args)
        params = self._get_joint_snv_mix_params(args)

        return JointSnvMixModel(priors, params, model=args.model)
    
    def _get_joint_snv_mix_priors(self, args):
        priors = JointSnvMixPriors()
        
        if args.mode == 'train' and args.priors_file is not None:
            priors.read_from_file(args.priors_file)
            
            print "Using custom initial priors from {0}".format(args.priors_file)
        
        return priors
    
    def _get_joint_snv_mix_params(self, args):
        params = JointSnvMixParameters()        
        
        if args.mode == 'train' and args.initial_parameters_file is not None:            
            params.read_from_file(args.initial_parameters_file)
            
            print "Using custom initial parameters from {0}".format(args.initial_parameters_file)
        
        if args.mode == 'classify' and args.parameters_file is not None:
            params.read_from_file(args.parameters_file)
            
            print "Using custom parameters from {0}".format(args.parameters_file)
        
        return params
