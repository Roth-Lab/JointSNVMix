from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import glob
import os

import Cython.Compiler.Options 
Cython.Compiler.Options.annotate = True

counter_includes = ['joint_snv_mix', 'joint_snv_mix/counters', 'include/samtools']

ref_iterator = Extension(
                         "joint_snv_mix.counters.ref_iterator",
                         ["joint_snv_mix/counters/ref_iterator.pyx"],
                         include_dirs=counter_includes
                         )
counter_row = Extension(
                        "joint_snv_mix.counters.counter_row",
                        ["joint_snv_mix/counters/counter_row.pyx"],
                        include_dirs=counter_includes
                        )

counter = Extension(
                    "joint_snv_mix.counters.counter",
                    ["joint_snv_mix/counters/counter.pyx"],
                    include_dirs=counter_includes
                    )

base_counter = Extension(
                        "joint_snv_mix.counters.base_counter",
                        ["joint_snv_mix/counters/base_counter.pyx"],
                        include_dirs=counter_includes
                        )
#
#quality_counter = Extension(
#                            "joint_snv_mix.counters.quality_counter",
#                            ["joint_snv_mix/counters/quality_counter.pyx"],
#                            include_dirs=counter_includes
#                            )
#
#joint_bin_counter = Extension(
#                            "joint_snv_mix.counters.joint_binary_counter",
#                            ["joint_snv_mix/counters/joint_binary_counter.pyx"],
#                            include_dirs=counter_includes
#                            )
#
#joint_quality_counter = Extension(
#                                  "joint_snv_mix.counters.joint_binary_quality_counter",
#                                  ["joint_snv_mix/counters/joint_binary_quality_counter.pyx"],
#                                  include_dirs=counter_includes
#                                  )
#
#positions_counter = Extension(
#                              "joint_snv_mix.counters.positions_counter",
#                              ["joint_snv_mix/counters/positions_counter.pyx"],
#                              include_dirs=counter_includes
#                              )
samtools_exclude = ("bamtk.c",
                    "razip.c",
                    "bgzip.c",
                    "main.c",
                    "calDepth.c",
                    "bam2bed.c",
                    "wgsim.c",
                    "md5fa.c",
                    "maq2sam.c")

all_samtools_files = ['joint_snv_mix/samtools.pyx']
all_samtools_files.extend(glob.glob('include/samtools/*.c'))
all_samtools_files.extend(glob.glob(os.path.join('include/samtools/', '*', '*.c')))

samtools_files = []

for file_name in all_samtools_files:
    if os.path.basename(file_name) in samtools_exclude:
        continue
    samtools_files.append(file_name)

samtools_includes = ['joint_snv_mix', 'include/samtools']

samtools = Extension(
                     "joint_snv_mix.samtools",
                     samtools_files,
                     include_dirs=samtools_includes,
                     libraries=[ "z", ]
                     )

utils_includes = [
                  'joint_snv_mix/utils'
                  ]

log_pdf = Extension(
                    "joint_snv_mix.utils.log_pdf",
                    ["joint_snv_mix/utils/log_pdf.pyx"],
                    include_dirs=utils_includes
                    )

ext_modules = [
               ref_iterator,
               counter,
               counter_row,
               base_counter,
               samtools,
               log_pdf
               ]

setup(
      name='JointSNVMix',
      version='0.8',
      description='A collection of tools for calling somatic mutations in paired tumour normal data.',
      author='Andrew Roth',
      author_email='andrewjlroth@gmail.com',
      url='http://compbio.bccrc.ca',
      
      packages=[ 
                'joint_snv_mix',
#                'joint_snv_mix.counters',
                'joint_snv_mix.utils'
                ],
      
      cmdclass={'build_ext': build_ext},
      ext_modules=ext_modules
     )
