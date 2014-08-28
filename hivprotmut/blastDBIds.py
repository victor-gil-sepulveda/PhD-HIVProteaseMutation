'''
Created on 28/8/2014

@author: victor
'''
import json
import sys
import hivprotmut.tools as tools
from hivprotmut.sequences.fastaFile import FastaFile
from hivprotmut.external.blast.blastpCommands import BlastpCommands

if __name__ == '__main__':
    parameters = json.loads(tools.remove_comments(open(sys.argv[1]).read()))
    
    # Prepare fasta file with input sequence
    input_fasta = FastaFile.open(parameters["query"]["fasta_file"])
    input_fasta.write("QUERY", parameters["query"]["sequence"])
    input_fasta.close()
    
    # Use blast
    alignments = BlastpCommands.find_closest_sequences(parameters["query"]["fasta_file"], 
                                                       parameters["blastp"])
    
    # Store ids
    handler = open(parameters["output"], "w")
    for alignment in alignments:
        handler.write(alignment["pdb"]["id"]+"\n")
    handler.close()
    
        