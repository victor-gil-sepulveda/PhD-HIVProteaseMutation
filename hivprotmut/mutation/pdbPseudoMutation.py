"""
Created on 14/8/2014

@author: victor
"""
import hivprotmut.tools as tools
import os

class PdbPseudoMutation(object):
    
    def __init__(self):
        pass
    
    @classmethod
    def get_mutated_residues(cls, master_sequence, target_sequence):
        """
        Precondition: SAME length
        TODO: generalize
        """
        mutations = {}
        for i in range(len(master_sequence)):
            if master_sequence[i] != target_sequence[i]:
                mutations[i+1] = tools.one_letter_residue_to_three(target_sequence[i]).upper()
        return mutations
    
    @classmethod
    def process_pdb(cls, mutations, input_pdb_path, output_pdb_path, parameters):
        handler = open(os.path.join("tmp_db",input_pdb_path),"r")
        mut_handler = open(os.path.join("tmp_db","%s.mut"%(output_pdb_path)),"w")
        cls.process_pdb_handlers(handler, mut_handler, mutations, parameters)
        handler.close()
        mut_handler.close()
    
    @classmethod   
    def process_pdb_handlers(cls, input_handler, output_handler, mutations, parameters):
        backbone_atoms = [' N  ', ' CA ', ' C  ', ' O  ', ' OXT']
        for line in input_handler:
            if line[:4] == "ATOM":
                atom = line[12:16]
                resnumber = int(line[23:26])
                if resnumber in mutations:
                    if atom in backbone_atoms:
                        output_handler.write(line[0:17]+mutations[resnumber]+line[20:])
                else:
                    output_handler.write(line)
            else:
                output_handler.write(line)
