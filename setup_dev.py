from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import glob
import os

import Cython.Compiler.Options 
Cython.Compiler.Options.annotate = True

#=======================================================================================================================
# Core
#=======================================================================================================================
core_includes = ['joint_snv_mix', 'joint_snv_mix/samtools', 'include/samtools']

core = []

for file_name in glob.glob("joint_snv_mix/*.pyx"):
    base_name = os.path.basename(file_name)
    root, ext = os.path.splitext(base_name)
    
    module = "joint_snv_mix.{0}".format(root)

    ext = Extension(module, [file_name, ], include_dirs=core_includes)
        
    core.append(ext)

#=======================================================================================================================
# Samtools
#=======================================================================================================================
samtools_files = ['bam.c',
                  'bam_aux.c',
                  'bam_import.c',
                  'bam_index.c',
                  'bam_pileup.c',
                  'bgzf.c',
                  'faidx.c',
                  'kstring.c',
                  'razf.c',
                  'sam.c',
                  'sam_header.c']
samtools_files = [os.path.join('include/samtools', x) for x in samtools_files]

samtools_includes = ['joint_snv_mix/samtools', 'include/samtools']

samtools = []

for file_name in glob.glob("joint_snv_mix/samtools/*.pyx"):
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

for file_name in glob.glob('joint_snv_mix/models/*.pyx'):
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
      
      cmdclass={'build_ext': build_ext},
      ext_modules=ext_modules
     )
