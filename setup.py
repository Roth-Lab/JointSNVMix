from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension( "joint_snv_mix.file_formats.pileup", ["joint_snv_mix/file_formats/pileup.pyx"] )]

setup( 
      name='JointSNVMix',
      version='0.5.0',
      description='Python SNV calling utility.',
      
      author='Andrew Roth',
      author_email='andrewjlroth@gmail.com',

      url='https://code.google.com/p/joint-snv-mix/',
      
      packages=[ 
                'joint_snv_mix',
                'joint_snv_mix.pre_processing',
                'joint_snv_mix.classification',
                'joint_snv_mix.post_processing',                
                ],
      
      cmdclass={'build_ext': build_ext},
      ext_modules=ext_modules,
      
      scripts=['jsm.py']
     )
