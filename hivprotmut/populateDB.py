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
from hivprotmut.external.blast.blastpCommands import BlastpCommands


if __name__ == '__main__':
    parameters = json.loads(tools.remove_comments(open(sys.argv[1]).read()))
    tools.create_folder(parameters["blast_database"]["structures"]["path"])
    fasta_db_filename = "%s.fas"%parameters["blast_database"]["blast"]["blast_database_name"]

    # Download pdbs that we do not have in our db    
    fasta_handler = FastaFile.open(fasta_db_filename,"w")
    ids_in_folder = tools.get_all_ids_from_folder(parameters["blast_database"]["structures"]["path"])
    for id in tools.get_all_ids_from_file(parameters["blast_database"]["structures"]["pdb_id_file"]):
        if id in ids_in_folder.keys():
            # The pdb is already in a database created USING THIS SCRIPT
            print "%s is already in the database."%id
            pdb = prody.parsePDB(os.path.join(parameters["blast_database"]["structures"]["path"], 
                                              ids_in_folder[id]))
        else:
            pdb, path = tools.get_pdb(id, parameters["blast_database"]["structures"]["download_selection"])
            os.system("mv %s %s"%(path,
                                  os.path.join(parameters["blast_database"]["structures"]["path"])))
        
        hw = prody.HierView(pdb.select("protein"))
        for chain in hw.iterChains():
            header = "%s:%s|PDBID|CHAIN|SEQUENCE"%(id, chain.getChid())
            fasta_handler.write(header, chain.getSequence())
            
    fasta_handler.close()
    
    # Create the initial blast database
    mask_db_filename = "%s.asnb"%parameters["blast_database"]["blast"]["blast_database_name"]
    BlastpCommands.create_database_from_fasta_file(fasta_db_filename,
                                                   mask_db_filename, 
                                                   parameters["blast_database"]["blast"])
    
    
    
    
        
        
        