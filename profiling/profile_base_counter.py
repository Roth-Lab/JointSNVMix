'''
Created on 2011-06-22

@author: Andrew Roth
'''
from argparse import Namespace

import pysam

from joint_snv_mix.pre_processing.bam_to_jcnt import paired_bams_to_jcnt

def run():
    args = Namespace()
    
    args.normal_bam_file_name = "/share/data/VBA0038/VBA0038.normal.bam"
    args.tumour_bam_file_name = "/share/data/VBA0038/VBA0038.tumour.bam"
    
    args.reference_genome_file_name = "/share/data/human_all.fa"
    
    args.jcnt_file_name = "/home/andrew/Desktop/test.jcnt"
    
    args.min_base_qual = 10
    args.min_map_qual = 10
    
    args.min_depth = 10
    
    args.positions_file = None
    
    paired_bams_to_jcnt(args)
            
if __name__ == "__main__":
#    cython_python()
#    visual_test()

#    import pstats, cProfile
#    
#    cProfile.runctx("cython_python()", globals(), locals(), "Profile.prof")
#    
#    s = pstats.Stats("Profile.prof")
#    s.strip_dirs().sort_stats("time").print_stats()
    
    import line_profiler
    
    from joint_snv_mix.pre_processing.bam_to_jcnt import write_all_positions
    
    profiler = line_profiler.LineProfiler(write_all_positions)
    
    run_str = "run()"
    
    profiler.run( run_str )
    profiler.print_stats()
    
#    import timeit
#    
#    python = "pure_python()"
#    cython = "cython_python()"
#    
#    t = timeit.Timer(python, "from __main__ import pure_python")
#    print t.timeit(10)
#    
#    t = timeit.Timer(cython, "from __main__ import cython_python")
#    print t.timeit(10)
#    run()
