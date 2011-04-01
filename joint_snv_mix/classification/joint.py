'''
Created on 2011-03-31

@author: Andrew Roth
'''
class JointModelRunner(ModelRunner):
    def _train(self, args):        
        if args.subsample_size > 0:
            counts = self._subsample(args.subsample_size)
        else:
            counts = self.reader.get_counts()
                   
        self.priors_parser.load_from_file(args.priors_file)
        self.priors = self.priors_parser.to_dict()
        
        self._write_priors()
        
        data = JointData(counts)
        
        self.parameters = self.model.train(data,
                                           self.priors,
                                           args.max_iters,
                                           args.convergence_threshold)

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
                              
            data = JointData(sub_counts)            
                
            resp = self.model.classify(data, self.parameters)
        
            self.writer.write_data(chr_name, sub_rows, resp)
            
            start = stop
            stop = min(stop + n, end)

class JointLatentVariables(EMLatentVariables):
    def _init_responsibilities(self, data):
        '''
        Intialise responsibilities via k-means clustering.
        '''
        a_1 = np.asarray(data.a['normal'], dtype=np.float64)
        b_1 = np.asarray(data.b['normal'], dtype=np.float64)
        p_1 = a_1 / (a_1 + b_1)
              
        a_2 = np.asarray(data.a['tumour'], dtype=np.float64)
        b_2 = np.asarray(data.b['tumour'], dtype=np.float64)
        p_2 = a_2 / (a_2 + b_2)

        shape = (data.nrows, 9)
        
        responsibilities = np.zeros(shape)
        
        init_centers = np.array((1., 0.5, 0.))
        
        cluster_centers_1, labels_1 = kmeans2(p_1, init_centers, minit='matrix')
        cluster_centers_2, labels_2 = kmeans2(p_2, init_centers, minit='matrix')

        labels = 3 * labels_1 + labels_2

        for id in range(9):
            index = labels == id
            
            responsibilities[index, id] = 1.
        
        self.responsibilities = responsibilities
        
class JointModelPriorParser(PriorParser):
    def __init__(self):
        PriorParser.__init__(self)
        
        self.ncomponent = self.nclass['normal'] * self.nclass['tumour']
    
    def _load_mix_weight_priors(self):       
        self.priors['kappa'] = np.zeros((self.ncomponent,))
            
        for i, genotype_tuple in enumerate(constants.joint_genotypes):
            genotype = "_".join(genotype_tuple)
        
            self.priors['kappa'][i] = self.parser.getfloat('kappa', genotype)
            
class JointParameterParser(ParameterParser):
    def __init__(self):
        ParameterParser.__init__(self)
        
        self.ncomponent = self.nclass['normal'] * self.nclass['tumour']

    def _load_mix_weights(self):       
        pi = np.zeros((self.ncomponent,))
            
        for i, genotype_tuple in enumerate(constants.joint_genotypes):
            genotype = "_".join(genotype_tuple)
        
            pi[i] = self.parser.getfloat('pi', genotype)
            
        self.parameters['pi'] = pi / pi.sum()            