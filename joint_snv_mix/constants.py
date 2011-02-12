'''
Created on 2010-05-24

@author: Andrew Roth
'''
#===============================================================================
# Genotypes
#===============================================================================
genotypes = ( "aa", "ab", "bb" )

joint_genotypes = []

for normal_genotype in genotypes:
    for tumour_genotype in genotypes:
        joint_genotypes.append( ( normal_genotype, tumour_genotype ) )

#===============================================================================
# Alphabet
#===============================================================================
nucleotides = ( 'A', 'C', 'G', 'T' )

#===============================================================================
# Genomes
#===============================================================================
genomes = ( "normal", "tumour" )

#=======================================================================================================================
# Multinomial Genotypes
#=======================================================================================================================
multinomial_genotypes = []

for i in range( 4 ):
    for j in range( i, 4 ):
        multinomial_genotypes.append( nucleotides[i] + nucleotides[j] )

joint_multinomial_genotypes = []

for normal_genotype in multinomial_genotypes:
    for tumour_genotype in multinomial_genotypes:
        joint_multinomial_genotypes.append( ( normal_genotype, tumour_genotype ) )

somatic_multinomial_genotypes_indices = []

for i, g in enumerate( joint_multinomial_genotypes ):
    # Check normal is homozygous.
    if g[0][0] != g[0][1]:
        continue
    
    if g[0] != g[1]:
        somatic_multinomial_genotypes_indices.append( i )
        
loh_multinomial_genotypes_indices = []

for i, g in enumerate( joint_multinomial_genotypes ):
    # Check normal is not homozygous.
    if g[0][0] == g[0][0]:
        continue
    
    # Check that tumour is homzygous
    if g[1][0] != g[1][1]:
        continue
    
    loh_multinomial_genotypes_indices.append( i )
    
matched_multinomial_genotypes_indices = []

for i, g in enumerate( joint_multinomial_genotypes ):
    # Check if genotypes match
    if g[0] != g[1]:
        continue
    
    matched_multinomial_genotypes_indices.append( i )
    
#=======================================================================================================================
# Extended multinomial
#=======================================================================================================================
extended_multinomial_genotypes = multinomial_genotypes[:]
extended_multinomial_genotypes.extend( ['ACG', 'ACT', 'AGT', 'CGT', 'ACGT'] ) 

joint_extended_multinomial_genotypes = []

for normal_genotype in extended_multinomial_genotypes:
    for tumour_genotype in extended_multinomial_genotypes:
        joint_extended_multinomial_genotypes.append( ( normal_genotype, tumour_genotype ) )


somatic_extended_multinomial_genotypes_indices = []

for i, g in enumerate( joint_extended_multinomial_genotypes ):
    # Check normal is diploid
    if len( g[0] ) != 2:
        continue
    
    # Check normal is homozygous.
    if g[0][0] != g[0][1]:
        continue
    
    if g[0] != g[1]:
        somatic_extended_multinomial_genotypes_indices.append( i )
        
loh_extended_multinomial_genotypes_indices = []

for i, g in enumerate( joint_extended_multinomial_genotypes ):
    # Check normal is diploid
    if len( g[0] ) != 2:
        continue
    
    # Check normal is not homozygous.
    if g[0][0] == g[0][0]:
        continue
    
    # Check that tumour is homzygous
    if g[1][0] != g[1][1]:
        continue
    
    loh_multinomial_genotypes_indices.append( i )
    
matched_extended_multinomial_genotypes_indices = []

for i, g in enumerate( joint_extended_multinomial_genotypes ):
    # Check if genotypes match
    if g[0] != g[1]:
        continue
    
    matched_multinomial_genotypes_indices.append( i )
    
#=======================================================================================================================
# Conan
#=======================================================================================================================
cn_state_map = {
                '1' : 3,
                '2' : 3,
                '3' : 3,
                '4' : 4,
                '5' : 5,
                '6' : 6
                }

conan_somatic_indices = {}

for cn_state in cn_state_map:
    nclass = cn_state_map[cn_state]
    
    conan_somatic_indices[cn_state] = []
    
    for i in range( 1, nclass ):
        conan_somatic_indices[cn_state].append( i )

binomial_alleles = ['a', 'b']

conan_genotypes = {}

for cn_state in cn_state_map:
    nclass = cn_state_map[cn_state]
    
    conan_genotypes[cn_state] = []
    
    for i in range( nclass + 1 ):
        genotype = 'a' * ( nclass - i ) + 'b' * i
        conan_genotypes[cn_state].append( genotype )
        
conan_joint_genotypes = {}        

for cn_state in cn_state_map:
    nclass = cn_state_map[cn_state]
    
    normal_genotypes = conan_genotypes['3']
    tumour_genotypes = conan_genotypes[cn_state]
    
    conan_joint_genotypes[cn_state] = []
    
    for ng in normal_genotypes:
        for tg in tumour_genotype:
            conan_joint_genotypes[cn_state].append( ( ng, tg ) )
        
    
    
