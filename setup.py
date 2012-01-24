from distutils.core import setup
from distutils.extension import Extension

import glob
import os

#=======================================================================================================================
# Counter
#=======================================================================================================================
counter_includes = ['joint_snv_mix', 'joint_snv_mix/samtools', 'include/samtools']

counter = Extension(
                    "joint_snv_mix.counter",
                    ["joint_snv_mix/counter.c"],
                    include_dirs=counter_includes
                    )

results = Extension(
                    "joint_snv_mix.results",
                    ["joint_snv_mix/results.c"],
                    include_dirs=counter_includes
                    )

run = Extension(
                "joint_snv_mix.run",
                ["joint_snv_mix/run.c"],
                include_dirs=counter_includes
                )

#=======================================================================================================================
# Samtools
#=======================================================================================================================
samtools_exclude = ("bamtk.c",
                    "razip.c",
                    "bgzip.c",
                    "main.c",
                    "calDepth.c",
                    "bam2bed.c",
                    "wgsim.c",
                    "md5fa.c",
                    "maq2sam.c")

all_samtools_files = glob.glob('include/samtools/*.c')
all_samtools_files.extend(glob.glob(os.path.join('include/samtools/', '*', '*.c')))

samtools_files = []

for file_name in all_samtools_files:
    if os.path.basename(file_name) in samtools_exclude:
        continue
    samtools_files.append(file_name)

samtools_includes = ['joint_snv_mix/samtools', 'include/samtools']

bam = Extension(
                 "joint_snv_mix.samtools.bam",
                 ['joint_snv_mix/samtools/bam.c', ] + samtools_files,
                 include_dirs=samtools_includes,
                 libraries=[ "z", ]
                 )

fasta = Extension(
                 "joint_snv_mix.samtools.fasta",
                 ['joint_snv_mix/samtools/fasta.c', ] + samtools_files,
                 include_dirs=samtools_includes,
                 libraries=[ "z", ]
                 )

pileup = Extension(
                   "joint_snv_mix.samtools.pileup",
                   ['joint_snv_mix/samtools/pileup.c', ] + samtools_files,
                   include_dirs=samtools_includes,
                   libraries=[ "z", ]
                   )

#=======================================================================================================================
# Models
#=======================================================================================================================
models_include = ['joint_snv_mix', 'joint_snv_mix/models', 'include/samtools']

models = Extension(
                    "joint_snv_mix.models.joint_snv_mix",
                    ["joint_snv_mix/models/joint_snv_mix.c"],
                    include_dirs=models_include
                    )

utils = Extension(
                  "joint_snv_mix.models.utils",
                  ["joint_snv_mix/models/utils.c"],
                  include_dirs=models_include
                  )


ext_modules = [
               bam,
               fasta,
               pileup,
               counter,
               models,
               utils,
               run,
               results               
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
                'joint_snv_mix.samtools',
                'joint_snv_mix.models'
                ],
      ext_modules=ext_modules,
      scripts=['jsm.py']
     )
