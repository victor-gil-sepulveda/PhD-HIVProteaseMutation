"""
Created on 13/8/2014

@author: victor
"""
from hivprotmut.sequences.fastaFile import FastaFile
import os
from hivprotmut.external.blast.blastOutputParser import BlastOutputParser
import hivprotmut.tools as tools

class BlastpCommands(object):
    """
    Namespace containing some blast related functions.
    """

    # This is an args template for a general blast query 
    BLASTP_ARGS =  " -query %s -out %s  -db %s -max_target_seqs %d -outfmt 5"

    def __init__(self):
        pass
    
    @classmethod
    def find_closest_sequences(cls, query_fasta_file, parameters):
        """
        Issues a blast command that queries a blast database in order to retrieve
        alignments with the closest sequences.
        
        :param query_fasta_file: the path of the fasta file with the query sequence.
        :param parameters: Dictionary holding the parameters.
        """
        command = "%s %s %s"%(parameters["exec"],
                              BlastpCommands.BLASTP_ARGS%(query_fasta_file, 
                                              parameters["blastp_output_file"], 
                                              parameters["search_db_name"], 
                                              parameters["max_target_sequences"]),
                              parameters["extra_args"])
        os.system(command)
        parser = BlastOutputParser(parameters["blastp_output_file"], True)
        return parser.alignments
    
    @classmethod
    def create_database_from_alignments(cls, alignments, parameters):
        """
        Function-like script to create a 
        """
        
        fasta_db_filename =  "%s.fas"%parameters["blast_database_name"]
        mask_db_filename = "%s.asnb"%parameters["blast_database_name"]
        
        fasta_db_handler = FastaFile.open(fasta_db_filename)
        for alignment in alignments:
            item_id = "%s_%s\n"%(alignment["pdb"]["id"], alignment["hit_chain"])
            fasta_db_handler.write(item_id, alignment["hit_seq__"])
        fasta_db_handler.close()
        
        cls.create_database_from_fasta_file(fasta_db_filename, mask_db_filename, parameters)
        
    @classmethod
    def create_database_from_fasta_file(cls, fasta_db_filename, mask_db_filename , parameters):       
        # Create new blast database
        os.system("%s -in %s -infmt fasta -outfmt maskinfo_asn1_bin -parse_seqids -out %s"%(
                                                parameters["segmasker_exe"],
                                                fasta_db_filename, 
                                                mask_db_filename))
        
        os.system("%s -in %s -dbtype prot -parse_seqids -mask_data %s -out %s -title %s"%(
                                       parameters["makeblastdb_exe"],
                                       fasta_db_filename,
                                       mask_db_filename,
                                       parameters["blast_database_name"],
                                       parameters["database_title"]))
        # Check
        os.system("%s -db %s -info"%(parameters["blastdbcmnd_exe"], 
                                     parameters["blast_database_name"])) 

        # Move everything to a folder (if requested)
        if parameters["blast_database_folder"] != "":
            tools.create_folder(parameters["blast_database_folder"])
            os.system("mv %s.* %s"%(
                                    parameters["blast_database_name"],
                                    parameters["blast_database_folder"])
                      )
