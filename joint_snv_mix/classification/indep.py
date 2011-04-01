'''
Created on 2011-03-31

@author: Andrew Roth
'''
class IndependentModelRunner(ModelRunner):
    def _train(self, args):
        if args.subsample_size > 0:
            counts = self._subsample(args.subsample_size)
        else:
            counts = self.reader.get_counts()
                   
        self.priors_parser.load_from_file(args.priors_file)
        self.priors = self.priors_parser.to_dict()
        
        self._write_priors()
        
        self.parameters = {}
        
        for genome in constants.genomes:
            data = IndependentData(counts, genome)
            
            self.parameters[genome] = self.model.train(data,
                                                       self.priors[genome],
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
            
            indep_resp = {}
            
            for genome in constants.genomes:                          
                data = IndependentData(sub_counts, genome)            
                
                indep_resp[genome] = self.model.classify(data, self.parameters[genome])
            
            joint_resp = self._get_joint_responsibilities(indep_resp)
        
            self.writer.write_data(chr_name, sub_rows, joint_resp)
            
            start = stop
            stop = min(stop + n, end)
            
    def _get_joint_responsibilities(self, resp):
        normal_resp = np.log(resp['normal'])
        tumour_resp = np.log(resp['tumour'])
        
        n = normal_resp.shape[0]
        
        nclass_normal = normal_resp.shape[1] 
        
        column_shape = (n, 1)
        
        log_resp = []
        
        for i in range(nclass_normal): 
            log_resp.append(normal_resp[:, i].reshape(column_shape) + tumour_resp)
        
        log_resp = np.hstack(log_resp)
        
        resp = np.exp(log_resp)
        
        return resp

class IndependentLatenVariables(EMLatentVariables):        
    def _init_responsibilities(self, data):
        a = np.asarray(data.a, dtype=np.float64)
        b = np.asarray(data.b, dtype=np.float64)
        
        shape = (data.nrows, 3)
        
        responsibilities = np.zeros(shape)
        
        p = a / (a + b)
        
        init_centers = np.array([1., 0.5, 0.])
        
        cluster_centers, labels = kmeans2(p, init_centers)
        
        sorted_centers = np.argsort(cluster_centers)
        
        for id in sorted_centers:
            index = labels == id
            
            responsibilities[index, id] = 1.0
        
        self.responsibilities = responsibilities

class IndepedendentPriorParser(PriorParser):    
    def _load_mix_weight_priors(self):       
        for genome in constants.genomes:            
            self.priors[genome]['kappa'] = np.zeros((self.nclass[genome],))
            
            for i, genotype in enumerate(constants.genotypes):                
                genome_genotype = "_".join((genome, genotype))
            
                self.priors[genome]['kappa'][i] = self.parser.getfloat('kappa', genome_genotype)
                
class IndepedendentParameterParser(ParameterParser):    
    def _load_mix_weights(self):       
        for genome in constants.genomes:            
            pi = np.zeros((self.nclass[genome],))
            
            for i, genotype in enumerate(constants.genotypes):                
                genome_genotype = "_".join((genome, genotype))
            
                pi[i] = self.parser.getfloat('pi', genome_genotype)
                
            self.parameters[genome]['pi'] = pi / pi.sum()