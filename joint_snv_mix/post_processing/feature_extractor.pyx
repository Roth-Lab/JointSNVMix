'''
Created on 2012-01-26

@author: Andrew Roth
'''
from collections import OrderedDict

cdef class FeatureExtractor(object):
    def __init__(self, BamFile normal_bam, BamFile tumour_bam):        
        self._normal_bam = normal_bam
        self._tumour_bam = tumour_bam
        
        priors = BinomialPriors()
        params = BinomialParameters()
        
        self._model = SnvMixTwoModel(priors, params)
    
    def get_features(self, JointBinaryCounterRow row):
        normal_features = self._get_genome_features(row, self._normal_bam)
        tumour_features = self._get_genome_features(row, self._tumour_bam)
        
        normal_features = tuple([('normal_' + x[0], x[1]) for x in normal_features])
        tumour_features = tuple([('tumour_' + x[0], x[1]) for x in tumour_features])
        
        ratio_features = self._get_ratio_features(OrderedDict(normal_features + tumour_features))
        
        jsm_features = self._get_joint_snv_mix_features(row)
        
        return OrderedDict(normal_features + tumour_features + ratio_features + jsm_features)
    
    cdef tuple _get_genome_features(self, JointBinaryCounterRow row, BamFile bam_file):
        cdef PileupIterator pileup_iter
    
        pileup_iter = bam_file.get_pileup_iterator(row._ref, start=row._pos, stop=row._pos + 1)
        
        while True:
            pileup_iter.advance_position()
            
            if pileup_iter._pos == row._pos:
                break
        
        pileup_column = pileup_iter.get_extended_pileup_column()
        
        features = self._parse_column(row, pileup_column)
        
        return features 
    
    cdef tuple _parse_column(self, JointBinaryCounterRow row, ExtendedPileupColumn pileup_column):
        site_features = self._get_site_features(pileup_column)
        
        ref_base = row._ref_base
        var_base = row._var_base
        
        ref_base_features = self._get_base_features(ref_base, pileup_column)
        var_base_features = self._get_base_features(var_base, pileup_column)        
        
        ref_base_features = tuple([('ref_' + x[0], x[1]) for x in ref_base_features])
        var_base_features = tuple([('var_' + x[0], x[1]) for x in var_base_features])
        
        return site_features + ref_base_features + var_base_features
    
    cdef tuple _get_joint_snv_mix_features(self, JointBinaryCounterRow row):
        snv_mix_probs = self._model.predict(row._data)
        
        somatic_prob_max = max(snv_mix_probs[1:3])
        non_somatic_prob_max = max(snv_mix_probs[0], max(snv_mix_probs[3:]))
        
        somatic_prob_sum = sum(snv_mix_probs[1:3])
        non_somatic_prob_sum = 1 - somatic_prob_sum
        
        wt_prob = snv_mix_probs[0]
        germline_prob = snv_mix_probs[4] + snv_mix_probs[8]
        loh_prob = snv_mix_probs[3] + snv_mix_probs[5]
        err_prob = snv_mix_probs[6] + snv_mix_probs[7]
        
        return (
                ('p_AA_AA', snv_mix_probs[0]),
                ('p_AA_AB', snv_mix_probs[1]),
                ('p_AA_BB', snv_mix_probs[2]),
                ('p_AB_AA', snv_mix_probs[3]),
                ('p_AB_AB', snv_mix_probs[4]),
                ('p_AB_BB', snv_mix_probs[5]),
                ('p_BB_AA', snv_mix_probs[6]),
                ('p_BB_AB', snv_mix_probs[7]),
                ('p_BB_BB', snv_mix_probs[8]),
                ('somatic_prob_max', somatic_prob_max),
                ('non_somatic_prob_max', non_somatic_prob_max),
                ('somatic_prob_sum', somatic_prob_sum),
                ('non_somatic_prob_sum', non_somatic_prob_sum),
                ('wt_prob', wt_prob),
                ('germline_prob', germline_prob),
                ('loh_prob', loh_prob),
                ('err_prob', err_prob)
                )
        
    #===================================================================================================================
    # Site specific features
    #===================================================================================================================
    cdef tuple _get_site_features(self, ExtendedPileupColumn pileup_column):
        depth = pileup_column._depth
        
        # Map qual features
        map_quals_rms = self._get_map_qual_rms(pileup_column)        
        map_quals_zeros = self._get_num_zero_map_qual_reads(pileup_column)
        
        # Alelle count
        alle_count_low_qual = self._get_total_allele_count(pileup_column, 0, 0)
        alle_count_high_qual = self._get_total_allele_count(pileup_column, 30, 30)
    
        return (
                ('depth', depth),
                ('map_quals_rms', map_quals_rms),
                ('map_quals_zeros', map_quals_zeros),
                ('alle_count_low_qual', alle_count_low_qual),
                ('alle_count_high_qual', alle_count_high_qual)
                )
        
    cdef double _get_map_qual_rms(self, ExtendedPileupColumn pileup_column):
        cdef int i
        cdef double map_qual_rms
        
        map_qual_rms = 0
    
        for i in range(pileup_column._depth):
            map_qual_rms += pow(pileup_column._map_quals[i], 2)
        
        return sqrt(map_qual_rms)
    
    cdef int _get_num_zero_map_qual_reads(self, ExtendedPileupColumn pileup_column):
        cdef int i, map_quals_zeros
        
        map_quals_zeros = 0
    
        for i in range(pileup_column._depth):
            if pileup_column._map_quals[i] == 0:
                map_quals_zeros += 1
        
        return map_quals_zeros
    
    cdef int _get_total_allele_count(self, ExtendedPileupColumn pileup_column, int base_qual, int map_qual):
        cdef char base[2], nucleotides[4]
        cdef int allele_count, i
        
        nucleotides[0] = 'A'
        nucleotides[1] = 'C'
        nucleotides[2] = 'G'
        nucleotides[3] = 'T'
        
        allele_count = 0
        
        for i in range(4):
            base[0] = nucleotides[i]
            base[1] = < char > NULL
            if pileup_column.get_nucleotide_count(base, base_qual, map_qual) > 0:
                allele_count += 1
        
        return allele_count
    
    #===================================================================================================================
    # Base specific features
    #===================================================================================================================
    cdef tuple _get_base_features(self, char * base, ExtendedPileupColumn pileup_column):
        # Allelic depth features
        base_count_low_quality = pileup_column.get_nucleotide_count(base, 0, 0)        
        base_count_high_quality = pileup_column.get_nucleotide_count(base, 30, 30)
    
        # Base quality features
        base_quals_sum = self._get_base_quals_sum(base, pileup_column, 1)                
        base_quals_square_sum = self._get_base_quals_sum(base, pileup_column, 2)
        
        # Mapping quality features
        map_quals_sum = self._get_map_quals_sum(base, pileup_column, 1)        
        map_quals_square_sum = self._get_map_quals_sum(base, pileup_column, 2)
        
        # Strand features
        forward_strand_count = self._get_forward_strand_count(base, pileup_column, 13)                        
        reverse_strand_count = pileup_column.get_nucleotide_count(base, 13, 0) - forward_strand_count 
        
        # Tail distrance features
        tail_distance_sum = self._get_tail_distance_sum(base, pileup_column, 1)        
        tail_distance_square_sum = self._get_tail_distance_sum(base, pileup_column, 2)
        
        # Homopolymer run before
        length_homopolymer_run_before_sum = self._get_base_matches_before_sum(base, pileup_column, 1)
        length_homopolymer_run_before_square_sum = self._get_base_matches_before_sum(base, pileup_column, 2)
        
        # Hompolymer run after
        length_homopolymer_run_after_sum = self._get_base_matches_after_sum(base, pileup_column, 1)
        length_homopolymer_run_after_square_sum = self._get_base_matches_after_sum(base, pileup_column, 2)
        
        # Length of homopolymer run containing base
        length_homopolymer_run_spanning_sum = self._get_base_matches_spanning_sum(base, pileup_column, 1)
        length_homopolymer_run_spanning_square_sum = self._get_base_matches_spanning_sum(base, pileup_column, 2)
                
        return (
                ('base_count_low_quality', base_count_low_quality),
                ('base_count_high_quality', base_count_high_quality),
                ('base_quals_sum', base_quals_sum),
                ('base_quals_square_sum', base_quals_square_sum),
                ('map_quals_sum', map_quals_sum),
                ('map_quals_square_sum', map_quals_square_sum),
                ('forward_strand_count', forward_strand_count),
                ('reverse_strand_count', reverse_strand_count),
                ('tail_distance_sum', tail_distance_sum),
                ('tail_distance_square_sum', tail_distance_square_sum),
                ('length_homopolymer_run_before_sum', length_homopolymer_run_before_sum),
                ('length_homopolymer_run_before_square_sum', length_homopolymer_run_before_square_sum),
                ('length_homopolymer_run_after_sum', length_homopolymer_run_after_sum),
                ('length_homopolymer_run_after_square_sum', length_homopolymer_run_after_square_sum),
                ('length_homopolymer_run_spanning_sum', length_homopolymer_run_spanning_sum),
                ('length_homopolymer_run_spanning_square_sum', length_homopolymer_run_spanning_square_sum)
                )
        
    cdef double _get_base_quals_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent):
        cdef int i
        cdef double base_qual_sum
        
        base_qual_sum = 0
    
        for i in range(pileup_column._depth):
            if base[0] == pileup_column._bases[i]:
                base_qual_sum += pow(pileup_column._base_quals[i], exponent)
        
        return base_qual_sum

    cdef double _get_map_quals_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent):
        cdef int i
        cdef double map_qual_sum
        
        map_qual_sum = 0
    
        for i in range(pileup_column._depth):
            if base[0] == pileup_column._bases[i]:
                map_qual_sum += pow(pileup_column._map_quals[i], exponent)
        
        return map_qual_sum
    
    cdef int _get_forward_strand_count(self, char * base, ExtendedPileupColumn pileup_column, int min_base_qual):
        cdef int i, forward_strand_count
        
        forward_strand_count = 0
        
        for i in range(pileup_column._depth):
            if (base[0] == pileup_column._bases[i]) and (pileup_column._base_quals[i] >= min_base_qual):
                if pileup_column._is_forward_strand[i]:
                    forward_strand_count += 1
        
        return forward_strand_count
    
    cdef double _get_tail_distance_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent):
        cdef int i
        cdef double tail_distance_sum
        
        tail_distance_sum = 0
        
        for i in range(pileup_column._depth):
            if base[0] == pileup_column._bases[i]:
                tail_distance_sum += pow(pileup_column._tail_distance[i], exponent)
        
        return tail_distance_sum
    
    cdef double _get_base_matches_before_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent):
        cdef int i
        cdef double run_length

        run_length = 0
        
        for i in range(pileup_column._depth):
            if base[0] == pileup_column._bases[i]:
                run_length += pow(pileup_column._base_matches_before[i], exponent)
        
        return run_length

    cdef double _get_base_matches_after_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent):
        cdef int ifeatures
        cdef double run_length

        run_length = 0
        
        for i in range(pileup_column._depth):
            if base[0] == pileup_column._bases[i]:
                run_length += pow(pileup_column._base_matches_after[i], exponent)
        
        return run_length        

    cdef double _get_base_matches_spanning_sum(self, char * base, ExtendedPileupColumn pileup_column, double exponent):
        cdef int i
        cdef double run_length

        run_length = 0
        
        for i in range(pileup_column._depth):
            if base[0] == pileup_column._bases[i]:
                run_length += pow(pileup_column._base_matches_before[i] + pileup_column._base_matches_after[i] + 1,
                                  exponent)
        
        return run_length
    
    #===================================================================================================================
    # Ratio Features
    #===================================================================================================================
    cdef tuple _get_ratio_features(self, features):
        cdef int normal_depth, tumour_depth
        
        normal_depth = features['normal_depth']
        tumour_depth = features['tumour_depth']
        
        # Forward strand ratios    
        ref_forward_strand_ratio = self._get_normalised_ratio(normal_depth,
                                                              tumour_depth,
                                                              features['normal_ref_forward_strand_count'],
                                                              features['tumour_ref_forward_strand_count'])
        
        var_forward_strand_ratio = self._get_normalised_ratio(normal_depth,
                                                              tumour_depth,
                                                              features['normal_var_forward_strand_count'],
                                                              features['tumour_var_forward_strand_count']) 

        # Reverse strand ratios
        ref_reverse_strand_ratio = self._get_normalised_ratio(normal_depth,
                                                              tumour_depth,
                                                              features['normal_ref_reverse_strand_count'],
                                                              features['tumour_ref_reverse_strand_count']) 

        var_reverse_strand_ratio = self._get_normalised_ratio(normal_depth,
                                                              tumour_depth,
                                                              features['normal_var_reverse_strand_count'],
                                                              features['tumour_var_reverse_strand_count']) 
        
        # Base qualities ratios
        ref_base_qualities_ratio = self._get_normalised_ratio(normal_depth,
                                                              tumour_depth,
                                                              features['normal_ref_base_quals_sum'],
                                                              features['tumour_ref_base_quals_sum'])      

        var_base_qualities_ratio = self._get_normalised_ratio(normal_depth,
                                                              tumour_depth,
                                                              features['normal_var_base_quals_sum'],
                                                              features['tumour_var_base_quals_sum'])          
        
        # Mapping qualities ratios
        ref_map_qualities_ratio = self._get_normalised_ratio(normal_depth,
                                                             tumour_depth,
                                                             features['normal_ref_map_quals_sum'],
                                                             features['tumour_ref_map_quals_sum'])      

        var_map_qualities_ratio = self._get_normalised_ratio(normal_depth,
                                                             tumour_depth,
                                                             features['normal_var_map_quals_sum'],
                                                             features['tumour_var_map_quals_sum'])

        # Tail_distance ratios
        ref_tail_distance_ratio = self._get_normalised_ratio(normal_depth,
                                                             tumour_depth,
                                                             features['normal_ref_tail_distance_sum'],
                                                             features['tumour_ref_tail_distance_sum'])      

        var_tail_distance_ratio = self._get_normalised_ratio(normal_depth,
                                                             tumour_depth,
                                                             features['normal_var_tail_distance_sum'],
                                                             features['tumour_var_tail_distance_sum'])
        
        # Tail_distance squared ratios
        ref_tail_distance_squared_ratio = self._get_normalised_ratio(normal_depth,
                                                             tumour_depth,
                                                             features['normal_ref_tail_distance_square_sum'],
                                                             features['tumour_ref_tail_distance_square_sum'])      

        var_tail_distance_squared_ratio = self._get_normalised_ratio(normal_depth,
                                                             tumour_depth,
                                                             features['normal_var_tail_distance_square_sum'],
                                                             features['tumour_var_tail_distance_square_sum'])
        
        # Base count low qual ratios
        ref_base_count_low_quality_ratio = self._get_normalised_ratio(normal_depth,
                                                                      tumour_depth,
                                                                      features['normal_ref_base_count_low_quality'],
                                                                      features['tumour_ref_base_count_low_quality'])      

        var_base_count_low_quality_ratio = self._get_normalised_ratio(normal_depth,
                                                                      tumour_depth,
                                                                      features['normal_var_base_count_low_quality'],
                                                                      features['tumour_var_base_count_low_quality'])

        # Base count high qual ratios
        ref_base_count_high_quality_ratio = self._get_normalised_ratio(normal_depth,
                                                                       tumour_depth,
                                                                       features['normal_ref_base_count_high_quality'],
                                                                       features['tumour_ref_base_count_high_quality'])      

        var_base_count_high_quality_ratio = self._get_normalised_ratio(normal_depth,
                                                                       tumour_depth,
                                                                       features['normal_var_base_count_high_quality'],
                                                                       features['tumour_var_base_count_high_quality'])
        
        # Low vs high qual ratios
        ref_normal_base_count_high_low_ratio = self._get_normalised_ratio(normal_depth,
                                                                          normal_depth,
                                                                          features['normal_ref_base_count_low_quality'],
                                                                          features['normal_ref_base_count_high_quality'])
        
        var_normal_base_count_high_low_ratio = self._get_normalised_ratio(normal_depth,
                                                                          normal_depth,
                                                                          features['normal_var_base_count_low_quality'],
                                                                          features['normal_var_base_count_high_quality'])

        ref_tumour_base_count_high_low_ratio = self._get_normalised_ratio(tumour_depth,
                                                                          tumour_depth,
                                                                          features['tumour_ref_base_count_low_quality'],
                                                                          features['tumour_ref_base_count_high_quality'])
        
        var_tumour_base_count_high_low_ratio = self._get_normalised_ratio(normal_depth,
                                                                          normal_depth,
                                                                          features['tumour_var_base_count_low_quality'],
                                                                          features['tumour_var_base_count_high_quality'])
        
        return (
                ('ref_forward_strand_ratio', ref_forward_strand_ratio),
                ('var_forward_strand_ratio', var_forward_strand_ratio),
                ('ref_reverse_strand_ratio', ref_reverse_strand_ratio),
                ('var_reverse_strand_ratio', var_reverse_strand_ratio),
                ('ref_base_qualities_ratio', ref_base_qualities_ratio),
                ('var_base_qualities_ratio', var_base_qualities_ratio),
                ('ref_map_qualities_ratio', ref_map_qualities_ratio),
                ('var_map_qualities_ratio', var_map_qualities_ratio),
                ('ref_tail_distance_ratio', ref_tail_distance_ratio),
                ('var_tail_distance_ratio', var_tail_distance_ratio),
                ('ref_tail_distance_squared_ratio', ref_tail_distance_squared_ratio),
                ('var_tail_distance_squared_ratio', var_tail_distance_squared_ratio),
                ('ref_base_count_low_quality_ratio', ref_base_count_low_quality_ratio),
                ('var_base_count_low_quality_ratio', var_base_count_low_quality_ratio),
                ('ref_base_count_high_quality_ratio', ref_base_count_high_quality_ratio),
                ('var_base_count_high_quality_ratio', var_base_count_high_quality_ratio),
                ('ref_normal_base_count_high_low_ratio', ref_normal_base_count_high_low_ratio),
                ('var_normal_base_count_high_low_ratio', var_normal_base_count_high_low_ratio),
                ('ref_tumour_base_count_high_low_ratio', ref_tumour_base_count_high_low_ratio),
                ('var_tumour_base_count_high_low_ratio', var_tumour_base_count_high_low_ratio)
                )
      
    
    cdef double _get_normalised_ratio(self,
                                      int normal_depth,
                                      int tumour_depth,
                                      double normal_feature,
                                      double tumour_feature):
        cdef double normal_normalised_feature, tumour_normalised_feature
        
        normal_normalised_feature = normal_feature / normal_depth
        tumour_normalised_feature = tumour_feature / tumour_depth
        
        if normal_normalised_feature == 0:
            if tumour_normalised_feature == 0:
                return 0
            else:
                return INFINITY
        else:        
            return tumour_normalised_feature / normal_normalised_feature 
        
