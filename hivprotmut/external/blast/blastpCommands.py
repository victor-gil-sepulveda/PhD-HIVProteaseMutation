"""
Created on 13/8/2014

@author: victor
"""
from hivprotmut.sequences.fastaFile import FastaFile
import os
from hivprotmut.external.blast.blastOutputParser import BlastOutputParser
from hivprotmut.tools import save_json

class BlastpCommands(object):

    def __init__(self):
        pass
    
    BLASTP_ARGS =  " -query %s -out %s  -db %s -max_target_seqs %d -outfmt 5"
    
    @classmethod
    def find_closest_sequences(cls, query_fasta_file, parameters):
        """
        """
        command = "%s %s %s"%(parameters["exec"],
                              BlastpCommands.BLASTP_ARGS%(query_fasta_file, 
                                              parameters["blastp_output_file"], 
                                              parameters["search_db_name"], 
                                              parameters["max_target_sequences"]),
                              parameters["extra_args"])
        print command
        os.system(command)
        parser = BlastOutputParser(parameters["blastp_output_file"])
        save_json(parser.alignments, parameters["alignments_file"])
        return parser.alignments
    
    
    @classmethod
    def create_database_from_alignments(cls, alignments, parameters):
        """
        """
        fasta_db_filename =  "%s.fasta"%parameters["new_database_name"]
        mask_db_filename = "%s.asnb"%parameters["new_database_name"]
        fasta_db_handler = FastaFile.open(fasta_db_filename)

        for alignment in alignments:
            item_id = "%s_%s\n"%(alignment["pdb"]["id"], alignment["hit_chain"])
            fasta_db_handler.write(item_id, alignment["hit_seq__"])
        fasta_db_handler.close()
            
        # Create new blast database
        os.system("segmasker -in %s -infmt fasta -outfmt maskinfo_asn1_bin -parse_seqids -out %s"%(fasta_db_filename, mask_db_filename))
        os.system("makeblastdb -in %s -dbtype prot -parse_seqids -mask_data %s -out %s -title %s"%(
                                       fasta_db_filename,
                                       mask_db_filename,
                                       parameters["new_database_name"],
                                       parameters["database_title"]))
        # Check
        os.system("blastdbcmd -db %s -info"%parameters["new_database_name"])