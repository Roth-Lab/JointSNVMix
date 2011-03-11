'''
Created on 2011-02-22

@author: Andrew Roth
'''
from joint_snv_mix.classification.data import MultinomialData
from joint_snv_mix.classification.model_runners import ModelRunner, ChromosomeModelRunner
from joint_snv_mix.classification.models import JointMultinomialModel
from joint_snv_mix.classification.parameter_parsers import JointMultinomialParameterParser
from joint_snv_mix.classification.prior_parsers import JointMultinomialPriorParser
from joint_snv_mix.file_formats.jmm import JointMultiMixWriter
from joint_snv_mix.file_formats.mcnt import MultinomialCountsReader

def run_multimix(args):
    if args.priors_file is None:
        args.train = False
    else:
        args.train = True
    
    if args.model == "joint":        
        runner = JointMultinomialRunner()    
    elif args.model == "chromosome":
        runner = ChromosomeMultinomialRunner()
            
    runner.run(args)

#=======================================================================================================================
# Runner
#=======================================================================================================================
class MultinomialModelRunner(ModelRunner):
    def run(self, args):
        self.reader = MultinomialCountsReader(args.mcnt_file_name)
        self.writer = JointMultiMixWriter(args.jmm_file_name)
        
        ModelRunner.run(self, args)
               
    def _train(self, args):
        if args.subsample_size > 0:
            counts = self._subsample(args.subsample_size)
        else:
            counts = self.reader.get_counts()
                   
        self.priors_parser.load_from_file(args.priors_file)
        self.priors = self.priors_parser.to_dict()
        
        self._write_priors()
        
        data = MultinomialData(counts)
        
        self.parameters = self.model.train(data, self.priors,
                                            args.max_iters, args.convergence_threshold)

    def _classify_chromosome(self, chr_name):
        counts = self.reader.get_counts(chr_name)
        jcnt_rows = self.reader.get_table(chr_name)
        
        end = self.reader.get_chr_size(chr_name)

        n = int(1e5)
        start = 0
        stop = min(n, end)
        

        while start < end:
            sub_counts = counts[start:stop]
            sub_rows = jcnt_rows[start:stop]
                              
            data = MultinomialData(sub_counts)            
                
            resp = self.model.classify(data, self.parameters)
        
            self.writer.write_data(chr_name, sub_rows, resp)
            
            start = stop
            stop = min(stop + n, end)
    
class JointMultinomialRunner(MultinomialModelRunner):
    def __init__(self):
        self.model = JointMultinomialModel()
        self.priors_parser = JointMultinomialPriorParser()
        self.parameter_parser = JointMultinomialParameterParser()
        
class ChromosomeMultinomialRunner(ChromosomeModelRunner):
    def __init__(self):
        self.data_class = MultinomialData
        
        self.model = JointMultinomialModel()
        self.priors_parser = JointMultinomialPriorParser()
        self.parameter_parser = JointMultinomialParameterParser()
        
    def run(self, args):
        self.reader = MultinomialCountsReader(args.jcnt_file_name)
        self.writer = JointMultiMixWriter(args.jsm_file_name)
        
        ModelRunner.run(self, args)
