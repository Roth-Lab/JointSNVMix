import joint_snv_mix.constants as constants

file_name = "joint_extended_multinomial_em.priors.cfg"

f = open( file_name, 'w' )

f.write( "[delta]" )
f.write( "\n" )

for mg in constants.multinomial_genotypes:
    if len( mg ) == 2:    
        for nuc in constants.nucleotides:
            if nuc in mg:
                f.write( "normal_{0}_{1} = {2}\n".format( mg, nuc, mg.count( nuc ) * 50 ) )
                f.write( "tumour_{0}_{1} = {2}\n".format( mg, nuc, mg.count( nuc ) * 50 ) )
                f.write( "\n" )
            else:
                f.write( "normal_{0}_{1} = 1\n".format( mg, nuc ) )
                f.write( "tumour_{0}_{1} = 1\n".format( mg, nuc ) )
            f.write( "\n" )
    elif len( mg ) == 3:
        for nuc in constants.nucleotides:
            if nuc in mg:
                f.write( "normal_{0}_{1} = {2}\n".format( mg, nuc, 1000 ) )
                f.write( "tumour_{0}_{1} = {2}\n".format( mg, nuc, 1000 ) )
                f.write( "\n" )
            else:
                f.write( "normal_{0}_{1} = 1\n".format( mg, nuc ) )
                f.write( "tumour_{0}_{1} = 1\n".format( mg, nuc ) )
            f.write( "\n" )
    elif len( mg ) == 4:
        for nuc in constants.nucleotides:
            if nuc in mg:
                f.write( "normal_{0}_{1} = {2}\n".format( mg, nuc, 1000 ) )
                f.write( "tumour_{0}_{1} = {2}\n".format( mg, nuc, 1000 ) )
                f.write( "\n" )
            else:
                f.write( "normal_{0}_{1} = 1\n".format( mg, nuc ) )
                f.write( "tumour_{0}_{1} = 1\n".format( mg, nuc ) )
            f.write( "\n" )


f.write( "\n" )
f.write( "[kappa]" )
f.write( "\n" )

for jmg in constants.joint_extended_multinomial_genotypes:
    jmg_string = "_".join( jmg )
    f.write( "{0} = 2\n".format( jmg_string ) )
    f.write( "\n" )
    
f.write( "scaling = 2\n" )
    
f.close()
