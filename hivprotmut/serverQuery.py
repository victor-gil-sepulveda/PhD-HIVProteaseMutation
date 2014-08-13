"""
Created on 13/8/2014

:author: victor

:brief: Given a query compares the sequences with the ones in our database and generates the mutant
"""
from hivprotmut.sequences.fastaFile import FastaFile
from hivprotmut.blast.blastpCommands import BlastpCommands

if __name__ == '__main__':
    seq = "PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
    
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
    BlastpCommands.find_closest_sequences("server_query.fasta", 
                                          blastp_parameters)
    
    # Get pdb from our database
    
    # 'Mutate' residues
    
    # Apply protein wizard
    
    
    
    
    