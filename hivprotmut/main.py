'''
Created on 12/8/2014

@author: victor
'''
import os
import prody
from hivprotmut.BlastOutputParser import BlastOutputParser

MAX_SEQUENCES = 50

###############################
### BLAST CONFIG
BLASTP_EXEC = "blastp"
BLASTP_ARGS =  " -query %s -out %s -remote -db %s -max_target_seqs %d -outfmt 5"
BLAST_DB = "pdbaa"
###############################


def find_closest_sequence_pdbs(fasta_file, output_file):
    """
    """
    os.system(BLASTP_EXEC + BLASTP_ARGS%(fasta_file, output_file, BLAST_DB, MAX_SEQUENCES))

def extract_pdbs_from_blast_output(blast_output_file):
    """
    """
    handler = open(blast_output_file,"r")
    pdb_ids = []
    
    for line in handler:
        id = line.split("|")[4]
        print id
        pdb_ids.append(id)
        
    handler.close()

if __name__ == '__main__':
    find_closest_sequence_pdbs("hivprotmut/test/HIV.fasta", "sequences.txt")
    extract_pdbs_from_blast_output("sequences.txt")
    bot =  BlastOutputParser("sequences.xml")
    bot.save("alignments.json")
    


