from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
               Extension("joint_snv_mix.pre_processing.base_counter", ["joint_snv_mix/pre_processing/base_counter.pyx"], include_dirs=['include/pysam', 'include/samtools'])
               ]

setup(
      name='JointSNVMix',
      version='0.6.3',
      description='Paired sample SNV calling utility.',
      
      author='Andrew Roth',
      author_email='andrewjlroth@gmail.com',

      url='https://code.google.com/p/joint-snv-mix/',
      
      packages=[ 
                'joint_snv_mix',
                'joint_snv_mix.classification',
                'joint_snv_mix.classification.utils',
                'joint_snv_mix.pre_processing',
                'joint_snv_mix.post_processing',
                'joint_snv_mix.file_formats'            
                ],
      
      scripts=['jsm.py'],
      cmdclass={'build_ext': build_ext},
      ext_modules=ext_modules
     )
