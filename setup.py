from distutils.core import setup

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
      
      scripts=['jsm.py']
     )
