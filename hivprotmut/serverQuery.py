"""
Created on 13/8/2014

:author: victor

:brief: Given a query compares the sequences with the ones in our database and generates the mutant
"""
import os
from hivprotmut.sequences.fastaFile import FastaFile
from hivprotmut.external.blast.blastpCommands import BlastpCommands
from hivprotmut.mutation.pdbPseudoMutation import PdbPseudoMutation
from hivprotmut.external.plop.plopCommands import PlopCommands
from hivprotmut import tools
import sys

if __name__ == '__main__':
    parameters = tools.remove_comments(open(sys.argv[1]).read())
    
    input_fasta = FastaFile.open(parameters["fasta_query"])
    input_fasta.write("QUERY", parameters["query"])
    input_fasta.close()
    
    alignments = BlastpCommands.find_closest_sequences(parameters["fasta_query"], 
                                          parameters["blastp"])
    
    # Get pdb from our database
    if len(alignments) >= 1:
        # Find mutated residues
        mutations = PdbPseudoMutation.get_mutated_residues(alignments[0]["query_seq"],
                                                    alignments[0]["hit_seq__"])

        # Get pdb path
        pdb_id = alignments[0]["pdb"]["id"].lower()
        processed_pdb_path = "%s.pdb.prot_lig_water"%pdb_id
        
        # "Mutate" pdb
        PdbPseudoMutation.process_pdb(mutations, 
                                      processed_pdb_path, 
                                      "mutated_tmp.pdb", 
                                      {})
        
        # Use plop to add the missing atoms
        plop_parameters = {
                           "exec": "/home/victor/plop/plop_ROT_tmp",
                           "control_file": "reconstruct_control_file"
        }
        PlopCommands.reconstruct_pdb("mutated_tmp.pdb", "mutated.pdb", plop_parameters)
        
        # Apply protein wizard to mutated pdb
        
        # Use prep wizard
        #$SCHRODINGER/utilities/prepwizard -noimpref -fix  to_mute.pdb muti.pdb
        
        # PDB is prepared for lig exploration!
    else:
        print "ERROR: impossible to find a suitable template from our database."
    
    
    
    
    
#         s = SequenceMatcher(None,
#                             alignments[0]["query_seq"],
#                             alignments[0]["hit_seq__"],
#                             autojunk=False)
#         print s.get_opcodes()

# Prody mutation
#         pdb = prody.parsePDB(os.path.join("tmp_db",processed_pdb_id))
#         heteros = pdb.select("hetero")
#         all_but_mutis = pdb.select("protein and not resnum "+" ".join([str(x) for x in mutations]))
#         backbone_mutis = pdb.select("backbone and resnum "+" ".join([str(x) for x in mutations]))
#         print all_but_mutis, backbone_mutis
#         prody.writePDB("kkk.pdb",(heteros+backbone_mutis+all_but_mutis))
#         