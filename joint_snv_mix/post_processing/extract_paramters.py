from joint_snv_mix.file_formats.jsm import JointSnvMixReader

def main( jsm_file_name ):
    reader = JointSnvMixReader( jsm_file_name )
    
    parameters = reader.get_parameters()
    
    for parameter_name, parameter_value in parameters.iteritems():
        if parameter_name == 'pi':
            continue
        
        print "Normal {0} : {1}".format( parameter_name, parameter_value[0] )
        print "Tumour {0} : {1}".format( parameter_name, parameter_value[1] )
    
    print "Mix-weights : {0}".format( parameters['pi'] )
    
    reader.close()
    

if __name__ == "__main__":
    import sys
    
    jsm_file_name = sys.argv[1]
    
    main( jsm_file_name )
