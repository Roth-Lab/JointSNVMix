import csv
import bisect
from collections import Counter

import pysam

from joint_snv_mix.file_formats.jcnt import JointCountsWriter
from joint_snv_mix.pre_processing.base_counter import BaseCounter

ascii_offset = 33

def bam_to_jcnt(args):
    if args.positions_file is not None:
        regions = convert_positions_to_regions(args.positions_file)
    else:
        regions = None
        
    normal_bam = pysam.Samfile(args.normal_bam_file_name, 'rb')
    tumour_bam = pysam.Samfile(args.tumour_bam_file_name, 'rb')
    
    normal_counter = BaseCounter(normal_bam, min_base_qual=args.min_base_qual, min_map_qual=args.min_base_qual)
    tumour_counter = BaseCounter(tumour_bam, min_base_qual=args.min_base_qual, min_map_qual=args.min_base_qual)
    
    ref_genome_fasta = pysam.Fastafile(args.reference_genome_file_name)
               
    converter = BamToJcntConverter(
                                   normal_counter,
                                   tumour_counter,
                                   ref_genome_fasta,
                                   min_depth=args.min_depth,
                                   regions=regions
                                   )
    
    converter.convert(args.jcnt_file_name)
                       
class BamToJcntConverter:
    def __init__(self, normal_counter, tumour_counter, ref_genome_fasta, min_depth=4, regions=None):
        self.normal_counter = normal_counter
        self.tumour_counter = tumour_counter
    
        self.ref_genome_fasta = ref_genome_fasta
        
        self.refs = sorted(set(self.normal_counter.refs) & set(self.tumour_counter.refs))
        
        self.min_depth = min_depth
        
        self.max_buffer = int(1e5)
        
        self.regions = regions
        
    def convert(self, jcnt_file_name):        
        with JointCountsWriter(jcnt_file_name) as writer:
            self.writer = writer
                                 
            if self.regions is None:
                self._convert_by_chrom()
            else:
                self._convert_by_region()
            
            self.writer.close()
        
    def _convert_by_chrom(self):
        for ref in self.refs:
            print ref
            
            joint_iter = JointPileupIterator(
                                       self.normal_counter.iter_ref(ref),
                                       self.tumour_counter.iter_ref(ref)
                                       )
            
            self._convert_iter(ref, joint_iter)
        
    def _convert_by_region(self):
        for ref in self.refs:
            if ref not in self.regions:
                continue
            
            print "Parsing {0} with {1} regions.".format(ref, len(self.regions[ref]))
            
            joint_iter = JointPileupRegionIterator(
                                                   self.normal_counter.iter_ref(ref),
                                                   self.tumour_counter.iter_ref(ref),
                                                   self.regions[ref]
                                                   )
            
            self._convert_iter(ref, joint_iter)
        
    def _convert_iter(self, ref, iter):        
        for normal_row, tumour_row in iter:
            if normal_row.depth < self.min_depth or tumour_row.depth < self.min_depth:
                continue
            
            pos = normal_row.position
            ref_base = self.ref_genome_fasta.fetch(ref, pos, pos + 1)        
            
            jcnt_entry = self._get_jcnt_entry(normal_row, tumour_row, pos, ref_base)

            if jcnt_entry is None:
                continue
            
            self.writer.add_row(ref, jcnt_entry)      
        
    def _get_jcnt_entry(self, normal_row, tumour_row, pos, ref_base):
        tumour_non_ref_base, tumour_counts = self._get_counts(ref_base, tumour_row.counts)        
        normal_non_ref_base, normal_counts = self._get_counts(ref_base, normal_row.counts, non_ref_base=tumour_non_ref_base)        
    
        # Check again for lines below read depth. The first check above speeds things up, though redundant.
        normal_depth = normal_counts[0] + normal_counts[1]
        tumour_depth = tumour_counts[0] + tumour_counts[1]
    
        if normal_depth < self.min_depth or tumour_depth < self.min_depth:
            return None
        
        jcnt_entry = [ pos, ref_base, normal_non_ref_base, tumour_non_ref_base ]
        jcnt_entry.extend(normal_counts)
        jcnt_entry.extend(tumour_counts)
        
        return jcnt_entry

    def _get_counts(self, ref_base, bases, non_ref_base=None):
        counter = Counter(bases)
        
        non_ref_base, counts = self._parse_counts(ref_base, counter, non_ref_base)
        
        return non_ref_base, counts
    
    def _parse_counts(self, ref_base, counter, non_ref_base=None):
        ref_counts = counter[ref_base]

        # Check if there is any non-ref bases.
        if non_ref_base is not None:
            non_ref_counts = counter[non_ref_base]
        else:            
            non_ref_base, non_ref_counts = counter.most_common(2)[1]
            
            if non_ref_counts == 0:
                non_ref_base = 'N'
        
        counts = (ref_counts, non_ref_counts)
        
        return non_ref_base, counts
    
class JointPileupIterator:
    def __init__(self, normal_iter, tumour_iter):
        self.normal_iter = normal_iter
        self.tumour_iter = tumour_iter
    
    def __iter__(self):
        return self
    
    def next(self):
        normal_counter_row = self.normal_iter.next()
        tumour_counter_row = self.tumour_iter.next()
        
        while True:                    
            normal_pos = normal_counter_row.position
            tumour_pos = tumour_counter_row.position
                  
            if normal_pos == tumour_pos:
                return normal_counter_row, tumour_counter_row
            elif normal_pos < tumour_pos:
                normal_counter_row = self.normal_iter.next()
            elif normal_pos > tumour_pos:
                tumour_counter_row = self.tumour_iter.next()
            else:
                raise Exception("Error in joint pileup iterator.")
            
class JointPileupRegionIterator:
    def __init__(self, normal_iter, tumour_iter, regions):
        self.normal_iter = normal_iter
        self.tumour_iter = tumour_iter
        
        self._init_index_lists(regions)
    
    def __iter__(self):
        return self
    
    def _init_index_lists(self, regions):
        starts = []
        ends = []
        
        for region in regions:
            start = region[0]
            end = region[1]
            
            starts.append(start)
            ends.append(end)
        

        self._starts = sorted(starts)
        self._ends = sorted(ends)
            
    def _in_regions(self, pos):
        # Before begining or after end of all regions.
        if pos < self._starts[0] or pos > self._ends[-1]:
            return False
        
        index = bisect.bisect_left(self._ends, pos)
                
        if self._starts[index] <= pos <= self._ends[index]:
            return True
        else:
            return False
    
    def next(self):
        normal_column = self.normal_iter.next()
        tumour_column = self.tumour_iter.next()
        
        while True:                    
            normal_pos = normal_column.pos
            tumour_pos = tumour_column.pos
                  
            if normal_pos == tumour_pos:
                if self._in_regions(normal_pos):                
                    return normal_column, tumour_column
                else:
                    normal_column = self.normal_iter.next()
                    tumour_column = self.tumour_iter.next()
            elif normal_pos < tumour_pos:
                normal_column = self.normal_iter.next()
            elif normal_pos > tumour_pos:
                tumour_column = self.tumour_iter.next()
            else:
                raise Exception("Error in joint pileup iterator.")

#=======================================================================================================================
# Functions
#=======================================================================================================================
def convert_positions_to_regions(positions_file_name):
    reader = csv.reader(open(positions_file_name), delimiter=' ')

    regions = {}
    
    row = reader.next()
    
    prev_chrom = row[0]
    prev_pos = int(row[1])
    
    region_chrom = prev_chrom
    region_start = prev_pos
    
    for row in reader:
        chrom = row[0]
        pos = int(row[1])
        
        if (chrom != prev_chrom) or (pos != prev_pos + 1):
            if region_chrom not in regions:
                regions[region_chrom] = []
            
            # Close current region and print
            region_end = prev_pos
            
            # Subtract 1 to make zero based.
            regions[region_chrom].append((region_start - 1, region_end - 1))
            
            # Open a new region
            region_chrom = chrom
            region_start = pos
        
        prev_chrom = chrom
        prev_pos = pos
    
    region_end = prev_pos
    regions[region_chrom].append((region_start - 1, region_end - 1))

    return regions
