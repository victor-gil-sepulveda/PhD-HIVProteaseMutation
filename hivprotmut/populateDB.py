'''
Created on 26/8/2014

@author: victor
'''
import hivprotmut.tools as tools
import sys
from hivprotmut.sequences.fastaFile import FastaFile
import json
import os
import prody

if __name__ == '__main__':
    parameters = json.loads(tools.remove_comments(open(sys.argv[1]).read()))
    
    tools.create_folder(parameters["blast_database"]["path"])
    
    #"HIV-1_Protease.PDBID.list")
    fasta_handler = FastaFile.open(parameters["blast_database"]["fasta"],"w")
    ids_in_folder = tools.get_all_ids_from_folder(parameters["blast_database"]["path"])
    
    for id in tools.get_all_ids_from_file(parameters["blast_database"]["pdb_id_file"]):
        if id in ids_in_folder.keys():
            # The pdb is already in a database created USING THIS SCRIPT
            pdb = prody.parsePDB(os.path.join(parameters["blast_database"]["path"], ids_in_folder[id]))
        else:
            pdb, _ = tools.get_pdb(id, parameters["blast_database"]["download_selection"])
        fasta_handler.write(id, tools.get_protein_sequence(pdb))
    fasta_handler.close()
    
    
        
        
        