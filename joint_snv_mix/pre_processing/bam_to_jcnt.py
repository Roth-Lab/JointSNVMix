import csv
from collections import Counter

import pysam

from joint_snv_mix.file_formats.jcnt import JointCountsWriter
import bisect

ascii_offset = 33

def bam_to_jcnt(args):
    if args.positions_file is not None:
        regions = convert_positions_to_regions(args.positions_file)
    else:
        regions = None
        
    normal_bam = pysam.Samfile(args.normal_bam_file_name, 'rb')
    tumour_bam = pysam.Samfile(args.tumour_bam_file_name, 'rb')
    
    ref_genome_fasta = pysam.Fastafile(args.reference_genome_file_name)
               
    converter = BamToJcntConverter(
                                   normal_bam,
                                   tumour_bam,
                                   ref_genome_fasta,
                                   min_depth=args.min_depth,
                                   min_bqual=args.min_base_qual,
                                   min_mqual=args.min_base_qual,
                                   regions=regions
                                   )
    
    converter.convert(args.jcnt_file_name)
                       
class BamToJcntConverter:
    def __init__(self, normal_bam, tumour_bam, ref_genome_fasta, min_depth=4, min_bqual=10, min_mqual=10, regions=None):
        self.normal_bam = normal_bam
        self.tumour_bam = tumour_bam
    
        self.ref_genome_fasta = ref_genome_fasta
        
        self.refs = sorted(set(self.normal_bam.references) & set(self.tumour_bam.references))
        
        self.min_depth = min_depth
        self.min_bqual = min_bqual
        self.min_mqual = min_mqual
        
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
                                       self.normal_bam.pileup(ref),
                                       self.tumour_bam.pileup(ref)
                                       )
            
            self._convert_iter(ref, joint_iter)
        
    def _convert_by_region(self):
        for ref in self.refs:
            if ref not in self.regions:
                continue
            
            print "Parsing {0} with {1} regions.".format(ref, len(self.regions[ref]))
            
            joint_iter = JointPileupRegionIterator(
                                                   self.normal_bam.pileup(ref),
                                                   self.tumour_bam.pileup(ref),
                                                   self.regions[ref]
                                                   )
            
            self._convert_iter(ref, joint_iter)
        
    def _convert_iter(self, ref, iter):        
        for normal_column, tumour_column in iter:
            if normal_column.n < self.min_depth or tumour_column.n < self.min_depth:
                continue
            
            pos = normal_column.pos
            ref_base = self.ref_genome_fasta.fetch(ref, pos, pos + 1)        
            
            jcnt_entry = self._get_jcnt_entry(normal_column, tumour_column, pos, ref_base)

            if jcnt_entry is None:
                continue
            
            self.writer.add_row(ref, jcnt_entry)      
        
    def _get_jcnt_entry(self, normal_column, tumour_column, pos, ref_base):
        normal_bases = self._parse_pileup_column(normal_column)
        tumour_bases = self._parse_pileup_column(tumour_column)
        
        tumour_non_ref_base, tumour_counts = self._get_counts(ref_base, tumour_bases)        
        normal_non_ref_base, normal_counts = self._get_counts(ref_base, normal_bases, non_ref_base=tumour_non_ref_base)        
    
        # Check again for lines below read depth. The first check above speeds things up, though redundant.
        normal_depth = normal_counts[0] + normal_counts[1]
        tumour_depth = tumour_counts[0] + tumour_counts[1]
    
        if normal_depth < self.min_depth or tumour_depth < self.min_depth:
            return None
    
        # Shift index to one based position.
        one_based_pos = pos + 1
        
        jcnt_entry = [ one_based_pos, ref_base, normal_non_ref_base, tumour_non_ref_base ]
        jcnt_entry.extend(normal_counts)
        jcnt_entry.extend(tumour_counts)
        
        return jcnt_entry
               
    def _parse_pileup_column(self, pileup_column):
        bases = []
        
        for read in pileup_column.pileups:
            if read.is_del:
                continue
                 
            mqual = read.alignment.mapq
            
            if mqual < self.min_mqual:
                continue            
            
            qpos = read.qpos
            
            bqual = ord(read.alignment.qual[qpos]) - ascii_offset            
            
            if bqual < self.min_bqual:
                continue
            
            base = read.alignment.seq[qpos]
            bases.append(base)
        
        return bases

    def _get_counts(self, ref_base, bases, non_ref_base=None):
        counter = Counter(bases)
        
        non_ref_base, counts = self._parse_counts(ref_base, counter, non_ref_base)
        
        return non_ref_base, counts
    
    def _parse_counts(self, ref_base, counter, non_ref_base=None):
        ref_counts = counter[ref_base]
        
        del counter[ref_base]
        del counter['N']
        
        # Check if there is any non-ref bases.
        if non_ref_base is not None:
            non_ref_counts = counter[non_ref_base]
        else:
            if len(counter) > 0:
                non_ref_base, non_ref_counts = counter.most_common(1)[0]
            else:
                non_ref_base = 'N'
                non_ref_counts = 0
        
        counts = (ref_counts, non_ref_counts)
        
        return non_ref_base, counts
    
class JointPileupIterator:
    def __init__(self, normal_iter, tumour_iter):
        self.normal_iter = normal_iter
        self.tumour_iter = tumour_iter
    
    def __iter__(self):
        return self
    
    def next(self):
        normal_column = self.normal_iter.next()
        tumour_column = self.tumour_iter.next()
        
        while True:                    
            normal_pos = normal_column.pos
            tumour_pos = tumour_column.pos
                  
            if normal_pos == tumour_pos:
                return normal_column, tumour_column
            elif normal_pos < tumour_pos:
                normal_column = self.normal_iter.next()
            elif normal_pos > tumour_pos:
                tumour_column = self.tumour_iter.next()
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
