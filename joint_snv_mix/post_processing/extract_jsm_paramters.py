from joint_snv_mix.file_formats.jsm import JointSnvMixReader

def extract_jsm_parameters(args):
    reader = JointSnvMixReader(args.jsm_file_name)    
    parameters = reader.get_parameters()    
    reader.close()
    
    recursive_print(parameters)

def recursive_print(params):
    for name, value in sorted(params.iteritems()):
        if isinstance(value, dict):
            print "-" * 100
            print name
            recursive_print(value)
        else:
            print name
            for row in value.tolist():
                print row

if __name__ == "__main__":
    import sys
    
    jsm_file_name = sys.argv[1]
    
    extract_jsm_parameters(jsm_file_name)
