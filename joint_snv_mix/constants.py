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

joint.append( 'junk' )

#===============================================================================
# Alphabet
#===============================================================================
nucleotides = ( 'A', 'C', 'G', 'T' )

#===============================================================================
# Genomes
#===============================================================================
genomes = ( "normal", "tumour", "junk" )

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
extended_multinomial_genotypes = multinomial_genotypes
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
    
