"""
Created on 13/8/2014

:author: victor

:brief: Given a query compares the sequences with the ones in our database and generates a "mutant"
structure using the closer structure as template.
"""

import sys
import json
from hivprotmut import tools
from hivprotmut.sequences.fastaFile import FastaFile
from hivprotmut.external.plop.plopCommands import PlopCommands
from hivprotmut.mutation.pdbPseudoMutation import PdbPseudoMutation
from hivprotmut.external.blast.blastpCommands import BlastpCommands
from hivprotmut.external.proteinwizard.pwCommands import PWCommands
import time
import glob
import os

if __name__ == '__main__':
    # Read params
    parameters = json.loads(tools.remove_comments(open(sys.argv[1]).read()))
    
    if "sequence" in parameters["query"]:
        # Prepare fasta file with input sequence
        input_fasta = FastaFile.open(parameters["query"]["fasta_file"])
        input_fasta.write("QUERY", parameters["query"]["sequence"])
        input_fasta.close()
    
    # Look for a match in our DB
    alignments = BlastpCommands.find_closest_sequences(parameters["query"]["fasta_file"], 
                                          parameters["blastp"])
    
    # Get pdb from our structure DB
    if len(alignments) >= 1:
        # Find mutated residues
        mutations = PdbPseudoMutation.get_mutated_residues(alignments[0]["query_seq"],
                                                    alignments[0]["hit_seq__"])

        # Get pdb path from structures DB
        pdb_id = alignments[0]["pdb"]["id"].lower()
        processed_pdb_path = os.path.join(parameters["query"]["ligand"],"%s.pdb"%pdb_id)
        
        # "Mutate" pdb
        PdbPseudoMutation.process_pdb(mutations, 
                                      processed_pdb_path, 
                                      parameters["mutation"])
        
        # Use PLOP to add the missing atoms
        PlopCommands.reconstruct_pdb(parameters["mutation"]["output_file"], 
                                     parameters["plop"])
        
        # Use prep wizard
        PWCommands.apply_protein_wizard(parameters["plop"]["output_file"], 
                                        parameters["protein_wizard"])

        # Monitor the creation of the file
        counter = 0
        while len(glob.glob('*.recover')) != 0:
            time.sleep(1)
            counter = counter + 1

        # PDB is prepared for lig exploration!
        print "[SUCCESS] This is the created pdb: %s"%(parameters["protein_wizard"]["output_file"])
        
    else:
        print "[ERROR] impossible to find a suitable template from our database."
