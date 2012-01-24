import os

root_dir = "../joint_snv_mix"

for dir_path, subFolders, files in os.walk(root_dir):
    for file_name in files:
        base_name, file_ext = os.path.splitext(file_name)
        
        file_name = os.path.join(dir_path, file_name)
        
        if file_ext in ['.c', '.so']:
            os.remove(file_name)
