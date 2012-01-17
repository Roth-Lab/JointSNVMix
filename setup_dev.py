from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import Cython.Compiler.Options 
Cython.Compiler.Options.annotate = True

#counter_includes = ['joint_snv_mix/counters', 'include/pysam', 'include/samtools']

#ref_iterator = Extension(
#                         "joint_snv_mix.counters.ref_iterator",
#                         ["joint_snv_mix/counters/ref_iterator.pyx"],
#                         include_dirs=counter_includes
#                         )
#counter_row = Extension(
#                        "joint_snv_mix.counters.counter_row",
#                        ["joint_snv_mix/counters/counter_row.pyx"],
#                        include_dirs=counter_includes
#                        )
#
#counter = Extension(
#                    "joint_snv_mix.counters.counter",
#                    ["joint_snv_mix/counters/counter.pyx"],
#                    include_dirs=counter_includes
#                    )
#
#base_counter = Extension(
#                        "joint_snv_mix.counters.base_counter",
#                        ["joint_snv_mix/counters/base_counter.pyx"],
#                        include_dirs=counter_includes
#                        )
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
samtools_files = ['joint_snv_mix/samtools.pyx', 'include/samtools/faidx.c']
samtools_includes = ['joint_snv_mix', 'include/samtools']

samtools = Extension(
                     "joint_snv_mix.samtools",
                     samtools_files,
                     include_dirs=samtools_includes
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
