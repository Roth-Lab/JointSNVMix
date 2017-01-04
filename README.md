JointSNVMix implements a probabilistic graphical model to analyse sequence data from tumour/normal pairs. 
The model draws statistical strength by analysing both genome jointly to more accurately classify germline and somatic mutations.

If you have questions see if they have been answered on the google group for JointSNVMix users at  http://groups.google.com/group/jointsnvmix-user-group. 
If you can't find an answer to your question, please post it on the google group.

For information about other software from the Shah lab see http://compbio.bccrc.ca/.

If you use JointSNVMix please cite :

[Andrew Roth et al., JointSNVMix : A Probabilistic Model For Accurate Detection Of Somatic Mutations In Normal/Tumour Paired Next Generation Sequencing Data](http://bioinformatics.oxfordjournals.org/content/early/2012/01/27/bioinformatics.bts053.abstract)

# Installing JointSNVMix

## Download

The source code implemented in Python is available under an open source license. 
The last version is available in dist/ folder of the git repository.

## Dependencies
Though not explicitly required by the JointSNVMix software, [samtools](http://www.htslib.org/) can be useful for creating bam, bam index and fasta index files that are required by JointSNVMix.

Before installing the software the following libraries are required.

  * [Python](http://www.python.org) >= 2.7

## Installing JointSNVMix
Once all dependencies have been installed JointSNVMix can be installed as follows.

  * tar -zxvf JointSNVMix-x.y.z.tar.gz
  * cd JointSNVMix-x.y.z
  * python setup.py install

There is also config/ folder alongside the setup.py file. 
This folder contains example priors, and parameters files for running the JointSNVMix software. 
It may be useful to copy this somewhere easily accessible.

# Running JointSNVMix

The JointSnvMix software package consists of a number of tools for calling somatic mutations in tumour/normal paired NGS data.
After installation the `jsm.py` command should be available on your system.
There are two subcommands `jsm.py classify` and `jsm.py train`. 
If you type `jsm.py classify -h` or `jsm.py train -h` the list of options along with help will be presented.


## Models
There are currently three models supported by both the `train` and `classify` commands. 
All models use the JointSNVMix mixture model which jointly analyses the normal and tumour genomes.
By default snvmix2 is used but other models can be specified using the `--model` flag.

## Description
  * binomial - Uses binomial densities in the mixture model this was previously referred to as the JointSnvMix1 mode.
  * snvmix2 - Uses snvmix2 densities in the mixture as described in the original SNVMix paper previously referred to as JointSnvMix2.
  * beta\_binomial - Uses beta-binomial densities in the mixture model new in version 0.8. The beta-binomial is a robust (in the statistical sense) alternative to binomial model. It can be beneficial when dealing with over-dispersed data. This is useful in cancer genomes since allelic frequencies at somatic mutations sites may deviate significantly from those expected under diploid model.
