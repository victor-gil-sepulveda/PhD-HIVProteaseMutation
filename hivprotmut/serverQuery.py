"""
Created on 13/8/2014

:author: victor

:brief: Given a query compares the sequences with the ones in our database and generates the mutant
"""
import os

###############################
### GLOBAL CONFIG
MAX_SEQUENCES = 1
###############################

###############################
### BLAST CONFIG
BLASTP_EXEC = "blastp"
BLASTP_ARGS =  " -query %s -out %s -remote -db %s -max_target_seqs %d -outfmt 5"
BLAST_DB = "hiv_filtered_sequences"
###############################

def find_closest_sequence_pdbs(fasta_file, output_file):
    """
    """
    os.system(BLASTP_EXEC + BLASTP_ARGS%(fasta_file, output_file, BLAST_DB, MAX_SEQUENCES))

if __name__ == '__main__':
    seq = "PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
    find_closest_sequence_pdbs()
    
    