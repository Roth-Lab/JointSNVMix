from distutils.core import setup
from distutils.extension import Extension

counter_includes = ['joint_snv_mix/counters', 'include/pysam', 'include/samtools']

ref_iterator = Extension(
                         "joint_snv_mix.counters.ref_iterator",
                         ["joint_snv_mix/counters/ref_iterator.c"],
                         include_dirs=counter_includes
                         )
counter_row = Extension(
                        "joint_snv_mix.counters.counter_row",
                        ["joint_snv_mix/counters/counter_row.c"],
                        include_dirs=counter_includes
                        )

counter = Extension(
                    "joint_snv_mix.counters.counter",
                    ["joint_snv_mix/counters/counter.c"],
                    include_dirs=counter_includes
                    )

base_counter = Extension(
                        "joint_snv_mix.counters.base_counter",
                        ["joint_snv_mix/counters/base_counter.c"],
                        include_dirs=counter_includes
                        )

quality_counter = Extension(
                            "joint_snv_mix.counters.quality_counter",
                            ["joint_snv_mix/counters/quality_counter.c"],
                            include_dirs=counter_includes
                            )

joint_bin_counter = Extension(
                            "joint_snv_mix.counters.joint_binary_counter",
                            ["joint_snv_mix/counters/joint_binary_counter.c"],
                            include_dirs=counter_includes
                            )

joint_quality_counter = Extension(
                                  "joint_snv_mix.counters.joint_binary_quality_counter",
                                  ["joint_snv_mix/counters/joint_binary_quality_counter.c"],
                                  include_dirs=counter_includes
                                  )

positions_counter = Extension(
                              "joint_snv_mix.counters.positions_counter",
                              ["joint_snv_mix/counters/positions_counter.c"],
                              include_dirs=counter_includes
                              )

classifier_includes = [
                       'joint_snv_mix/counters',
                       'joint_snv_mix/utils',
                       'joint_snv_mix/classifiers'
                       ]
classifier_includes.extend(counter_includes)

classifier = Extension(
                        "joint_snv_mix.classifiers.classifier",
                        ["joint_snv_mix/classifiers/classifier.c"],
                        include_dirs=classifier_includes
                        )

indep_fisher_classifier = Extension(
                                    "joint_snv_mix.classifiers.independent_fisher",
                                    ["joint_snv_mix/classifiers/independent_fisher.c"],
                                    include_dirs=classifier_includes
                                    )

joint_fisher_classifier = Extension(
                                    "joint_snv_mix.classifiers.joint_fisher",
                                    ["joint_snv_mix/classifiers/joint_fisher.c"],
                                    include_dirs=classifier_includes
                                    )

threshold_classifier = Extension(
                                 "joint_snv_mix.classifiers.threshold",
                                 ["joint_snv_mix/classifiers/threshold.c"],
                                 include_dirs=classifier_includes
                                 )

snv_mix_classifier = Extension(
                                 "joint_snv_mix.classifiers.snv_mix",
                                 ["joint_snv_mix/classifiers/snv_mix.c"],
                                 include_dirs=classifier_includes
                                 )

snv_mix_qualities_classifier = Extension(
                                         "joint_snv_mix.classifiers.snv_mix_qualities",
                                         ["joint_snv_mix/classifiers/snv_mix_qualities.c"],
                                         include_dirs=classifier_includes
                                         )

joint_snv_mix_classifier = Extension(
                                     "joint_snv_mix.classifiers.joint_snv_mix",
                                     ["joint_snv_mix/classifiers/joint_snv_mix.c"],
                                     include_dirs=classifier_includes
                                     )

joint_snv_mix_qualities_classifier = Extension(
                                               "joint_snv_mix.classifiers.joint_snv_mix_qualities",
                                               ["joint_snv_mix/classifiers/joint_snv_mix_qualities.c"],
                                               include_dirs=classifier_includes
                                               )


classifiers_shared = Extension(
                              "joint_snv_mix.classifiers.shared",
                              ["joint_snv_mix/classifiers/shared.c"],
                              include_dirs=classifier_includes
                              )

utils_includes = [
                  'joint_snv_mix/utils'
                  ]

fisher_exact_test = Extension(
                                "joint_snv_mix.utils.fisher_exact_test",
                                ["joint_snv_mix/utils/fisher_exact_test.c"],
                                include_dirs=utils_includes
                                )

special_functions = Extension(
                              "joint_snv_mix.utils.special_functions",
                              ["joint_snv_mix/utils/special_functions.c"],
                              include_dirs=utils_includes
                              )

log_pdf = Extension(
                    "joint_snv_mix.utils.log_pdf",
                    ["joint_snv_mix/utils/log_pdf.c"],
                    include_dirs=utils_includes
                    )

trainers_include = [
                    'joint_snv_mix/trainers',
                    ]
trainers_include.extend(counter_includes)

snv_mix_trainer = Extension(
                            "joint_snv_mix.trainers.snv_mix",
                            ["joint_snv_mix/trainers/snv_mix.c"],
                            include_dirs=trainers_include
                            )

joint_snv_mix_trainer = Extension(
                                  "joint_snv_mix.trainers.joint_snv_mix",
                                  ["joint_snv_mix/trainers/joint_snv_mix.c"],
                                  include_dirs=trainers_include
                                  )
ext_modules = [
               counter,
               ref_iterator,
               counter_row,
               base_counter,
               quality_counter,
               joint_bin_counter,
               joint_quality_counter,
               classifiers_shared,
               classifier,
               indep_fisher_classifier,
               joint_fisher_classifier,
               threshold_classifier,
               snv_mix_classifier,
               snv_mix_qualities_classifier,
               joint_snv_mix_classifier,
               joint_snv_mix_qualities_classifier,
               log_pdf,
               snv_mix_trainer,
               joint_snv_mix_trainer,
               positions_counter,
               fisher_exact_test,
               special_functions
               ]

setup(
      name='JointSNVMix',
      version='0.7.2',
      description='A collection of tools for calling somatic mutations in paired tumour normal data.',
      author='Andrew Roth',
      author_email='andrewjlroth@gmail.com',
      url='http://compbio.bccrc.ca',
      
      packages=[ 
                'joint_snv_mix',
                'joint_snv_mix.classifiers',
                'joint_snv_mix.counters',
                'joint_snv_mix.runners',
                'joint_snv_mix.trainers',
                'joint_snv_mix.utils'
                ],
      ext_modules=ext_modules,
      scripts=['jsm.py']
     )
