"""
Created on 13/8/2014

:author: victor

:brief: Given a query compares the sequences with the ones in our database and generates the mutant
"""
from hivprotmut.sequences.fastaFile import FastaFile
from hivprotmut.external.blast.blastpCommands import BlastpCommands
import os
from difflib import SequenceMatcher

if __name__ == '__main__':
    #        X            X                                                        X
    seq = "PQLTLWQRPLVTIKICGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAILTVLVGPTPVNIIGRNLLTQIGCTLNF"
    
    input_fasta = FastaFile.open("server_query.fasta")
    input_fasta.write("TEST", seq)
    input_fasta.close()
    
    blastp_parameters = {
                        "exec": "blastp",
                        "extra_args": "",
                        "blastp_output_file": "squery_sequence.xml",
                        "search_db_name": "newdb",
                        "max_target_sequences": 1,
                        "alignments_file": "squery_alignments.json"
                        }
    alignments = BlastpCommands.find_closest_sequences("server_query.fasta", 
                                          blastp_parameters)
    
    # Get pdb from our database
    if len(alignments) >= 1:
        # Find mutated residues
        mutations = []
        for i in range(len(alignments[0]["query_seq"])):
            if alignments[0]["query_seq"][i] != alignments[0]["hit_seq__"][i]:
                mutations.append(i+1)
#         s = SequenceMatcher(None,
#                             alignments[0]["query_seq"],
#                             alignments[0]["hit_seq__"],
#                             autojunk=False)
#         print s.get_opcodes()
        
        # Leave backbone only of mutated residues
        pdb_id = alignments[0]["pdb"]["id"].lower()
        processed_pdb_id = "%s.pdb.prot_lig_water"%pdb_id
        handler = open(os.path.join("tmp_db",processed_pdb_id),"r")
        mut_handler = open(os.path.join("tmp_db","%s.mut"%(processed_pdb_id)),"w")
        backbone_atoms = [' N  ', ' CA ', ' C  ', ' O  ', ' OXT']
        for line in handler:
            if line[:4] == "ATOM":
                atom = line[12:16]
                resnumber = int(line[23:26])
                if atom in backbone_atoms:
                    mut_handler.write(line)
                else:
                    if not resnumber in mutations:
                        mut_handler.write(line)
        handler.close()
        mut_handler.close()
        # Apply protein wizard
        
        os.system("$")
        
        # Use plop to add the missing atoms
        
        # Use prep wizard
        #$SCHRODINGER/utilities/prepwizard -noimpref -fix  to_mute.pdb muti.pdb
        
        # PDB is prepared for lig exploration!
    else:
        print "ERROR: impossible to find a suitable template from our database."
    
    
    
    
    