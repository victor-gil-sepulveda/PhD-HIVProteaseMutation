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
    tools.create_folder(parameters["structures"]["path"])
    fasta_db_filename = "%s.fas"%parameters["blast"]["blast_database_name"]

    # Download pdbs that we do not have in our db    
    fasta_handler = FastaFile.open(fasta_db_filename,"w")
    for pdb_id in tools.get_ids_from_file_list(parameters["structures"]["pdb_id_files"]):
        pdb, path = tools.get_pdb_from_remote_or_db(pdb_id, 
                                                    parameters["structures"]["download_selection"],
                                                    parameters["structures"]["path"])
        
        hw = prody.HierView(pdb.select("protein"))
        for chain in hw.iterChains():
            header = "%s:%s|PDBID|CHAIN|SEQUENCE"%(pdb_id, chain.getChid())
            fasta_handler.write(header, chain.getSequence())
        
        # Get rid of the pdb
        os.system("rm %s"%path)
            
    fasta_handler.close()
    
    # Create the initial blast database
    mask_db_filename = "%s.asnb"%parameters["blast"]["blast_database_name"]
    BlastpCommands.create_database_from_fasta_file(fasta_db_filename,
                                                   mask_db_filename, 
                                                   parameters["blast"])
    
    
    
    
        
        
        