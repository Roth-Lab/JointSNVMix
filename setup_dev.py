from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import Cython.Compiler.Options 
Cython.Compiler.Options.annotate = True

counter_includes = ['joint_snv_mix/counters', 'include/pysam', 'include/samtools']

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

joint_bin_counter = Extension(
                            "joint_snv_mix.counters.joint_binary_counter",
                            ["joint_snv_mix/counters/joint_binary_counter.pyx"],
                            include_dirs=counter_includes
                            )

classifier_includes = [
                       'joint_snv_mix/counters',
                       'joint_snv_mix/utils',
                       'joint_snv_mix/classifiers'
                       ]
classifier_includes.extend(counter_includes)

base_classifier = Extension(
                            "joint_snv_mix.classifiers.classifier",
                            ["joint_snv_mix/classifiers/classifier.pyx"],
                            include_dirs=classifier_includes
                            )

indep_fisher_classifier = Extension(
                                    "joint_snv_mix.classifiers.independent_fisher",
                                    ["joint_snv_mix/classifiers/independent_fisher.pyx"],
                                    include_dirs=classifier_includes
                                    )

joint_fisher_classifier = Extension(
                                    "joint_snv_mix.classifiers.joint_fisher",
                                    ["joint_snv_mix/classifiers/joint_fisher.pyx"],
                                    include_dirs=classifier_includes
                                    )

threshold_classifier = Extension(
                                 "joint_snv_mix.classifiers.threshold",
                                 ["joint_snv_mix/classifiers/threshold.pyx"],
                                 include_dirs=classifier_includes
                                 )

snv_mix_classifier = Extension(
                                 "joint_snv_mix.classifiers.snv_mix",
                                 ["joint_snv_mix/classifiers/snv_mix.pyx"],
                                 include_dirs=classifier_includes
                                 )

joint_snv_mix_classifier = Extension(
                                     "joint_snv_mix.classifiers.joint_snv_mix",
                                     ["joint_snv_mix/classifiers/joint_snv_mix.pyx"],
                                     include_dirs=classifier_includes
                                     )

utils_includes = [
                  'joint_snv_mix/utils'
                  ]

fisher_exact_test = Extension(
                                "joint_snv_mix.utils.fisher_exact_test",
                                ["joint_snv_mix/utils/fisher_exact_test.pyx"],
                                include_dirs=classifier_includes
                                )

special_functions = Extension(
                              "joint_snv_mix.utils.special_functions",
                              ["joint_snv_mix/utils/special_functions.pyx"],
                              include_dirs=classifier_includes
                              )

ext_modules = [ 
               counter,
               base_counter,
               joint_bin_counter,
               base_classifier,
               indep_fisher_classifier,
               joint_fisher_classifier,
               threshold_classifier,
               snv_mix_classifier,
               joint_snv_mix_classifier,
               fisher_exact_test,
               special_functions
               ]

setup(
      name='JointSNVMix',
      version='0.7',
      description='A collection of tools for calling somatic mutations in paired tumour normal data.',
      author='Andrew Roth',
      author_email='andrewjlroth@gmail.com',
      url='http://compbio.bccrc.ca',
      
      packages=[ 
                'bam_counter'              
                ],
      
      cmdclass={'build_ext': build_ext},
      ext_modules=ext_modules
     )
