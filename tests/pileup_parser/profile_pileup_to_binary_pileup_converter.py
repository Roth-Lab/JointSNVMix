'''
Created on 2010-07-30

@author: Andrew Roth
'''
from pyleup.pileup.binary_pileup import PileupToBinaryPileupConverter, BinaryPileupWriter

import line_profiler

input_file = "../../../data/ex1.pileup"
output_file = "../../../data/profile.bpil"

converter = PileupToBinaryPileupConverter()

run_str = "converter.convert_file( input_file, output_file )"

profiler = line_profiler.LineProfiler( PileupToBinaryPileupConverter.convert_file,
                                       BinaryPileupWriter.write_columns )

profiler.run( run_str )
profiler.print_stats()

import pstats, cProfile

cProfile.runctx(run_str, globals(), locals(), "Profile.prof")

s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()