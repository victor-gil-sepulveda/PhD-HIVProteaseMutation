'''
Created on 27/8/2014

@author: victor
'''
import sys
import os
import hivprotmut.tools as tools

if __name__ == '__main__':
    mandatory_list = sys.argv[1]
    origin_folder = sys.argv[2]
    send_to_folder = sys.argv[3]
    list_ids = tools.get_ids_from_file(mandatory_list)
    for pdb_id in list_ids:
        in_path = os.path.join(origin_folder,"%s.pdb"%pdb_id)
        out_path = os.path.join(send_to_folder,"%s.pdb"%pdb_id)
        os.system("cp %s %s"%(in_path, out_path))
