"""
Created on 1/9/2014

@author: victor
"""
import sys
import os
import glob
import prody 
from hivprotmut.structures.pdbcuration import CurationSelections

if __name__ == '__main__':
    final_db_folder = sys.argv[1]
    com_file = sys.argv[2]
    
    com_handler = open(com_file, "w")
    ligand_folders = os.listdir(final_db_folder) # first level are ligands
    txt_root = os.path.split(final_db_folder)[1]
    for path in ligand_folders:
        files = glob.glob(os.path.join(final_db_folder, path, "*.pdb"))
        for pdb_file in files:
            pdb = prody.parsePDB(pdb_file)
            txt_pdb = os.path.split(pdb_file)[1]
            ligand = pdb.select(CurationSelections.HEAVY_LIGAND_SELECTION)  
            com = prody.calcCenter(ligand)
            txt_pdb_file = os.path.join(path, txt_pdb)
            com_handler.write("%s %.3f %.3f %.3f\n"%(
                                                   txt_pdb_file,
                                                   com[0],
                                                   com[1],
                                                   com[2]                                                   
                                                   ))
    com_handler.close()