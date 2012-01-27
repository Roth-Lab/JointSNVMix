from distutils.core import setup
from distutils.extension import Extension

import glob
import os

#=======================================================================================================================
# Core
#=======================================================================================================================
core_includes = ['joint_snv_mix', 'joint_snv_mix/samtools', 'include/samtools']

core = []

for file_name in glob.glob("joint_snv_mix/*.c"):
    base_name = os.path.basename(file_name)
    root, ext = os.path.splitext(base_name)
    
    module = "joint_snv_mix.{0}".format(root)
    
    ext = Extension(module, [file_name, ], include_dirs=core_includes)
        
    core.append(ext)

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

samtools = []

for file_name in glob.glob("joint_snv_mix/samtools/*.c"):
    base_name = os.path.basename(file_name)
    root, ext = os.path.splitext(base_name)
    
    module = "joint_snv_mix.samtools.{0}".format(root)
    
    ext = Extension(module, [file_name, ] + samtools_files, include_dirs=samtools_includes, libraries=[ "z", ])
        
    samtools.append(ext)

#=======================================================================================================================
# Models
#=======================================================================================================================
models_include = ['joint_snv_mix', 'joint_snv_mix/models', 'include/samtools']

models = []

for file_name in glob.glob('joint_snv_mix/models/*.c'):
    base_name = os.path.basename(file_name)
    root, ext = os.path.splitext(base_name)
    
    module = "joint_snv_mix.models.{0}".format(root)
    
    ext = Extension(module, [file_name, ], include_dirs=models_include)
    
    models.append(ext)

ext_modules = []
ext_modules.extend(core)
ext_modules.extend(models)
ext_modules.extend(samtools)

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
