"""
Created on 12/8/2014

@author: victor
"""
import os
from hivprotmut.BlastOutputParser import BlastOutputParser
import hivprotmut.tools as tools
import shutil

###############################
### GLOBAL CONFIG
MAX_SEQUENCES = 50
###############################

###############################
### BLAST CONFIG
BLASTP_EXEC = "blastp"
BLASTP_ARGS =  " -query %s -out %s -remote -db %s -max_target_seqs %d -outfmt 5"
BLAST_DB = "pdbaa"
###############################

###############################
### DB CONFIG
PDB_TMP_DATABASE_FOLDER = "tmp_db"
PDB_DATABASE_FOLDER = "db"
###############################

###############################
### PREPROCESSING CONFIG
SELECTION_STRING = "protein"
###############################


def find_closest_sequence_pdbs(fasta_file, output_file):
    """
    """
    os.system(BLASTP_EXEC + BLASTP_ARGS%(fasta_file, output_file, BLAST_DB, MAX_SEQUENCES))

if __name__ == '__main__':
    tools.create_folder(PDB_TMP_DATABASE_FOLDER)
    tools.create_folder(PDB_DATABASE_FOLDER)
    
    # Find close sequences 
#     find_closest_sequence_pdbs("hivprotmut/test/HIV.fasta", "sequences.xml")
    
    # Parse results
    bot =  BlastOutputParser("sequences.xml")
    bot.save("alignments.json")
    print "Found %d alignments"%(len(bot.alignments))

    # Get the ids of pdbs without gap (backbones must be equal)
    filtered_alignments = [alignment_info for alignment_info in bot.alignments if alignment_info["gaps"] == 0]
    # IMPROVEMENT : LEAVE STRUCTURES WHERE THE GAPS ARE AT THE BEGGINING OR THE END TO SOME EXTENT
    print "We have filtered %d structures because the alignment had gaps."%(len(bot.alignments) - len(filtered_alignments))
    
    # Get the pdbs
    for alignment in filtered_alignments:
        pdb, pdb_path = tools.get_pdb(alignment["hit_id"], SELECTION_STRING)
        tmp_path = os.path.join(PDB_TMP_DATABASE_FOLDER, os.path.basename(pdb_path))
        shutil.move(pdb_path, tmp_path)
        
        
        
        


