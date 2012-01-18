'''
Created on 2011-08-11

@author: Andrew Roth
'''
import csv
    
#=======================================================================================================================
# Functions for running classification.
#=======================================================================================================================
def classify(args):
    parser_factory = ParserFactory()
    classifier_factory = ClassifierFactory()
    
    if args.model == "joint_snv_mix_one":
        parser = parser_factory.get_joint_snv_mix_one_parser(args)
        classifier = classifier_factory.get_joint_snv_mix_two_classifier(args)
    elif args.model == "joint_snv_mix_two":
        parser = parser_factory.get_joint_snv_mix_one_parser(args)
        classifier = classifier_factory.get_joint_snv_mix_two_classifier(args)
    else:
        raise Exception("{0} is not a valid model.".format(args.model))
    
    classify_data_set(parser, classifier, args)

def classify_data_set(parser, classifier, args):
    if args.out_file == '-':
        writer = StdoutWriter()
    else:
        writer = FileWriter(args.out_file)
    
    if args.chrom is not None:
        positions_iter = parser.get_chrom_iterator(args.chrom)
    else:
        positions_iter = parser.get_genome_iterator()
    
    for position in positions_iter:
        if not position.has_variant:
            continue 
        
        probs = classifier.predict_probabilities(position.data)
        
        writer.write_position(position, probs)

class FileWriter(object):
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
    
    def __init__(self, file_name):
        self._writer = csv.DictWriter(open(file_name, 'w'), FileWriter.fields, delimiter='\t')
        
        self._writer.writeheader()
        
    def write_position(self, position, probs):
        row = {}
        row['chrom'] = position.ref
        row['position'] = position.coord
        row['ref_base'] = position.ref_base
        row['var_base'] = position.var_base
        
        counts = position.counts
        
        row['normal_counts_a'] = counts[0]
        row['normal_counts_b'] = counts[1]
        row['tumour_counts_a'] = counts[2]
        row['tumour_counts_b'] = counts[3]
        
        probs = dict(zip(FileWriter.probs, probs))
        
        row = dict(row.items() + probs.items())
        
        self._writer.write(row)

class StdoutWriter(object):
    def write_position(self, position, probs):
        out_row = [
                   position.ref,
                   position.coord,
                   position.ref_base,
                   position.var_base                   
                   ]
        
        out_row.extend(position.counts)
        out_row.extend(probs)
        
        print "\t".join(out_row) 
         
#=======================================================================================================================
# Functions for training a model.
#=======================================================================================================================
def train(args):
    parser_factory = ParserFactory()
    classifier_factory = ClassifierFactory()
    
    if args.model == "joint_snv_mix_one":
        parser = parser_factory.get_joint_snv_mix_one_parser(args)
        classifier = classifier_factory.get_joint_snv_mix_two_classifier(args)
    elif args.model == "joint_snv_mix_two":
        parser = parser_factory.get_joint_snv_mix_one_parser(args)
        classifier = classifier_factory.get_joint_snv_mix_two_classifier(args)
    else:
        raise Exception("{0} is not a valid model.".format(args.model))
    
    data_set = create_training_data_set(parser, args)
    
    classifier.fit(data_set)
    
    classifier.params.save_to_file(args.params_file)        

def create_training_data_set(parser, args):
    i = 0
    
    data_set = []
    
    if args.chrom is not None:
        positions_iter = parser.get_chrom_iterator(args.chrom)
    else:
        positions_iter = parser.get_genome_iterator()
    
    for position in positions_iter:
        if args.min_depth <= position.depth <= args.max_depth:
            i += 1
            
            if i % args.skip_size == 0:
                data_set.append(position.data)
    
    return data_set

#=======================================================================================================================
# Helper factory classes for setting up parsers and classifiers.
#=======================================================================================================================
class ParserFactory(object):
    def get_joint_snv_mix_one_parser(self, args):
        parser = JointSnvMixOneDataParser(args.normal_file,
                                          args.tumour_file,
                                          args.min_base_qual,
                                          args.min_map_qual)
        
        return parser
    
    def get_joint_snv_mix_two_parser(self, args):
        parser = JointSnvMixTwoDataParser(args.normal_file,
                                          args.tumour_file)
        
        return parser

class ClassifierFactory(object):
    def get_joint_snv_mix_one_classifier(self, args):
        priors = self._get_joint_snv_mix_priors(args)
        params = self._get_joint_snv_mix_params(args)

        return JointSnvMixOneClassifier(priors=priors, params=params)
    
    def get_joint_snv_mix_two_classifier(self, args):
        priors = self._get_joint_snv_mix_priors(args)
        params = self._get_joint_snv_mix_params(args)

        return JointSnvMixTwoClassifier(priors=priors, params=params)
    
    def _get_joint_snv_mix_priors(self, args):
        priors = JointSnvMixPriors()
        
        if args.priors_file is not None:
            priors.load_from_file(args.priors_file)
            
            print "Using custom initial priors from {0}".format(args.init_params_file)
        else:
            print "Using default priors"
        
        print priors
        
        return priors
    
    def _get_joint_snv_mix_params(self, args):
        params = JointSnvMixParameters()        
        
        if args.init_params_file is not None:            
            params.load_from_file(args.init_params_file)
            
            print "Using custom initial parameters from {0}".format(args.init_params_file)
        else:
            print "Using default initial parameters."
        
        print params
        
        return params
