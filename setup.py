from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension( "joint_snv_mix.utils.pileup_parser", ["lib/joint_snv_mix/utils/pileup_parser.pyx"] )]

setup( 
      name='JointSNVMix',
      version='0.4.0',
      description='Python pileup parsing and SNV calling utilities.',
      
      author='Andrew Roth',
      author_email='andrewjlroth@gmail.com',

      url='http://gitorious.org/pyleup/joint_snv_mix',
      
      package_dir = {'': 'lib'},
      
      
      packages=[ 
                'joint_snv_mix',
                'joint_snv_mix.file_formats',
                'joint_snv_mix.snv_tools',
                'joint_snv_mix.snv_tools.joint_model',
                'joint_snv_mix.utils'                
                ],
      
      cmdclass={'build_ext': build_ext},
      ext_modules=ext_modules,
      
      scripts=[
               'bin/run_joint_snv_mix.py',
               'bin/call_joint_genotypes.py',
               'bin/pileup_to_jcnt.py',
               'bin/extract_jsm_positions.py'
               ]
     )
